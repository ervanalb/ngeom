use enum_dispatch::enum_dispatch;
use ngeom::ops::*;
use ngeom::re2::{AntiEven, AntiScalar, Vector};
use ngeom::scalar::*;
use ngeom_polygon::graph::EdgeSet;
use std::ops::{Add, Index as StdIndex, Mul};

pub trait Space {
    type Scalar: Copy
        + Ring
        + Rational
        + Sqrt<Output = Self::Scalar>
        + Trig<Output = Self::Scalar>
        + From<f32>;
    type AntiScalar: Copy;
    type Vector: Copy
        + Transform<Self::AntiEven, Output = Self::Vector>
        + Mul<Self::Scalar, Output = Self::Vector>
        + Add<Output = Self::Vector>;
    type AxisOfRotation: Copy
        + Mul<Self::Scalar, Output = Self::AxisOfRotation>
        + AntiMul<Self::AntiScalar, Output = Self::AxisOfRotation>;
    type AntiEven: Copy + AxisAngle<Self::AxisOfRotation, Self::Scalar>;
}

pub trait Index: Copy + Ord + Eq {}

impl Index for usize {}

#[derive(Clone, Copy, Debug)]
pub enum Dir {
    Fwd,
    Rev,
}

pub trait VertexCollection: StdIndex<Self::Index, Output = <Self::Space as Space>::Vector> {
    type Space: Space;
    type Index: Index;
}

pub trait PointStream {
    type Point;

    fn push(&mut self, point: Self::Point);
}

impl<POINT> PointStream for Vec<POINT> {
    type Point = POINT;

    fn push(&mut self, point: POINT) {
        self.push(point);
    }
}

#[derive(Clone)]
pub struct FaceBoundary<EI: Index> {
    pub edges: Vec<(EI, Dir)>,
}

#[derive(Clone)]
pub struct Face<EI: Index> {
    pub boundaries: Vec<FaceBoundary<EI>>,
}

pub trait EdgeCollection:
    StdIndex<Self::Index, Output = Edge<Self::Space, Self::VertexIndex>>
{
    type Space: Space;
    type Index: Index;
    type VertexIndex: Index;

    fn iter(&self) -> impl Iterator<Item = &Edge<Self::Space, Self::VertexIndex>>;
    //fn iter_outgoing(
    //    &self,
    //    vi: Self::VertexIndex,
    //) -> impl Iterator<Item = &Edge<Self::Space, Self::VertexIndex>>;
}

pub trait FaceCollection {
    type Space: Space;
    type Index: Index;
    type EdgeIndex: Index;
    type VertexIndex: Index;

    fn iter(&self) -> impl Iterator<Item = &Face<Self::EdgeIndex>>;
}

#[enum_dispatch]
pub trait EdgeTrait<SPACE: Space, VI: Index> {
    fn endpoints(&self) -> Option<(VI, VI)>;
    fn x(
        &self,
        vertices: &impl VertexCollection<Space = SPACE, Index = VI>,
        t: SPACE::Scalar,
    ) -> SPACE::Vector;
    fn interpolate(
        &self,
        vertices: &impl VertexCollection<Space = SPACE, Index = VI>,
        dir: Dir,
        point_stream: &mut impl PointStream<Point = SPACE::Vector>,
    );
}

#[derive(Clone, Debug)]
pub struct Line<VI: Index> {
    pub start: VI,
    pub end: VI,
}

impl<SPACE: Space, VI: Index> EdgeTrait<SPACE, VI> for Line<VI> {
    fn endpoints(&self) -> Option<(VI, VI)> {
        Some((self.start, self.end))
    }

    fn x(
        &self,
        vertices: &impl VertexCollection<Space = SPACE, Index = VI>,
        t: SPACE::Scalar,
    ) -> SPACE::Vector {
        vertices[self.start] * (SPACE::Scalar::one() - t) + vertices[self.end] * t
    }

    fn interpolate(
        &self,
        vertices: &impl VertexCollection<Space = SPACE, Index = VI>,
        dir: Dir,
        point_stream: &mut impl PointStream<Point = SPACE::Vector>,
    ) {
        point_stream.push(
            vertices[match dir {
                Dir::Fwd => self.start,
                Dir::Rev => self.end,
            }],
        );
    }
}

#[derive(Clone, Debug)]
pub struct Arc<SPACE: Space, VI: Index> {
    pub start: VI,
    pub axis: SPACE::AxisOfRotation,
    pub end_angle: SPACE::Scalar,
    pub end: VI,
}

impl<SPACE: Space, VI: Index> EdgeTrait<SPACE, VI> for Arc<SPACE, VI> {
    fn endpoints(&self) -> Option<(VI, VI)> {
        Some((self.start, self.end))
    }

    fn x(
        &self,
        vertices: &impl VertexCollection<Space = SPACE, Index = VI>,
        t: SPACE::Scalar,
    ) -> SPACE::Vector {
        let start_pt = vertices[self.start];
        start_pt.transform(SPACE::AntiEven::axis_angle(self.axis, self.end_angle * t))
    }

    fn interpolate(
        &self,
        vertices: &impl VertexCollection<Space = SPACE, Index = VI>,
        dir: Dir,
        point_stream: &mut impl PointStream<Point = SPACE::Vector>,
    ) {
        for i in 0..10 {
            // TODO dynamic spacing
            let i = match dir {
                Dir::Fwd => i,
                Dir::Rev => 9 - i,
            };
            let t = SPACE::Scalar::from(i as f32 / 10.);
            point_stream.push(self.x(vertices, t))
        }
    }
}

//#[derive(Clone, Debug)]
//pub struct Circle<SPACE: Space> {
//    pub axis: SPACE::AxisOfRotation,
//    pub pt: SPACE::Vector,
//}

#[enum_dispatch(EdgeTrait<SPACE, VI>)]
#[derive(Clone)]
pub enum Edge<SPACE: Space, VI: Index> {
    Line(Line<VI>),
    Arc(Arc<SPACE, VI>),
    //Circle(Circle<SPACE>),
    //CubicBezier-- or NURBS?
    //OffsetCubicBezier-- or NURBS?
}

#[derive(Debug)]
pub enum PatchTopologyError {
    Empty,
    Open,
    Branching,
}

pub struct Interpolation<POINT, UV> {
    pub points: Vec<POINT>,
    pub uv: Vec<UV>,
    pub edges: EdgeSet,
}

pub fn interpolate<
    UV,
    VC: VertexCollection,
    EC: EdgeCollection<Space = VC::Space, VertexIndex = VC::Index>,
    FC: FaceCollection<Space = VC::Space, EdgeIndex = EC::Index, VertexIndex = VC::Index>,
    IntoUvFn: Fn(<VC::Space as Space>::Vector) -> UV,
>(
    vertices: &VC,
    edges: &EC,
    faces: &FC,
    into_uv: IntoUvFn,
) -> Interpolation<<VC::Space as Space>::Vector, UV> {
    // First, interpolate the boundaries of the faces
    // and store their connectivity
    let mut points = Vec::<<VC::Space as Space>::Vector>::new();
    let mut polyedges = EdgeSet::new();

    for face in faces.iter() {
        for boundary in face.boundaries.iter() {
            let start = points.len();
            for &(edge_index, dir) in boundary.edges.iter() {
                edges[edge_index].interpolate(vertices, dir, &mut points);
            }
            let end = points.len();
            polyedges.insert_loop(start..end);
        }
    }

    // Now, re-interpret those boundary points into (U, V) coordinates.
    let uv: Vec<_> = points.iter().cloned().map(into_uv).collect();

    Interpolation {
        points,
        uv,
        edges: polyedges,
    }
}

pub struct Space2D;
impl Space for Space2D {
    type Scalar = f32;
    type AntiScalar = AntiScalar<f32>;
    type Vector = Vector<f32>;
    type AxisOfRotation = Vector<f32>;
    type AntiEven = AntiEven<f32>;
}

#[derive(Clone)]
pub struct VecVertex(pub Vec<Vector<f32>>);

impl StdIndex<usize> for VecVertex {
    type Output = Vector<f32>;

    fn index(&self, idx: usize) -> &Vector<f32> {
        self.0.index(idx)
    }
}

impl VertexCollection for VecVertex {
    type Space = Space2D;
    type Index = usize;
}

pub struct VecEdge(pub Vec<Edge<Space2D, usize>>);

impl StdIndex<usize> for VecEdge {
    type Output = Edge<Space2D, usize>;

    fn index(&self, idx: usize) -> &Edge<Space2D, usize> {
        self.0.index(idx)
    }
}

impl EdgeCollection for VecEdge {
    type Space = Space2D;
    type Index = usize;
    type VertexIndex = usize;

    fn iter(&self) -> impl Iterator<Item = &Edge<Space2D, usize>> {
        self.0.iter()
    }
    //fn iter_outgoing(&self, _vi: usize) -> impl Iterator<Item = &Edge<Space2D, usize>> {
    //    self.0.iter().filter(move |e| e.start == vi)
    //}
}

pub struct VecFace(pub Vec<Face<usize>>);

impl FaceCollection for VecFace {
    type Space = Space2D;
    type Index = usize;
    type EdgeIndex = usize;
    type VertexIndex = usize;

    fn iter(&self) -> impl Iterator<Item = &Face<usize>> {
        self.0.iter()
    }
    //fn iter_outgoing(&self, _vi: usize) -> impl Iterator<Item = &Edge<Space2D, usize>> {
    //    self.0.iter().filter(move |e| e.start == vi)
    //}
}
