use enum_dispatch::enum_dispatch;
use ngeom::ops::*;
use ngeom::scalar::*;
use std::collections::BTreeSet;
use std::fmt::Debug;
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

pub trait EdgeCollection {
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
trait EdgeTrait<SPACE: Space, VI: Index> {
    fn endpoints(&self) -> Option<(VI, VI)>;
    fn x(
        &self,
        vertices: &impl VertexCollection<Space = SPACE, Index = VI>,
        t: SPACE::Scalar,
    ) -> SPACE::Vector;
    fn interpolate(
        &self,
        vertices: &impl VertexCollection<Space = SPACE, Index = VI>,
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
        point_stream: &mut impl PointStream<Point = SPACE::Vector>,
    ) {
        point_stream.push(vertices[self.start]);
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
        point_stream: &mut impl PointStream<Point = SPACE::Vector>,
    ) {
        for i in 0..10 {
            // TODO dynamic spacing
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

#[enum_dispatch(EdgeTrait)]
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

pub fn triangulate<
    VC: VertexCollection,
    EC: EdgeCollection<Space = VC::Space, VertexIndex = VC::Index>,
    FC: FaceCollection<Space = VC::Space, EdgeIndex = EC::Index, VertexIndex = VC::Index>,
>(
    vertices: &VC,
    edges: &EC,
    faces: &FC,
) -> Result<Vec<<VC::Space as Space>::Vector>, PatchTopologyError> {
    let start_vi = edges.iter().next().ok_or(PatchTopologyError::Empty)?.start;
    let mut cur_vi = start_vi;
    let mut visited = BTreeSet::<VC::Index>::new();
    let mut polygon = Vec::<<VC::Space as Space>::Vector>::new();
    loop {
        // Find the next edge to interpolate
        let mut possible_next_vis = edges.iter_outgoing(cur_vi);
        let cur_edge = possible_next_vis.next().ok_or(PatchTopologyError::Open)?;
        if possible_next_vis.next().is_some() {
            // More than one edge returned
            return Err(PatchTopologyError::Branching);
        }
        if visited.contains(&cur_edge.end) {
            return Err(PatchTopologyError::Branching);
        }

        cur_edge.interpolate(&vertices, &mut polygon);
        // TODO
        //println!(
        //    "Interpolate edge! {:?} from {:?} to {:?}",
        //    cur_edge.curve, vertices[cur_edge.start], vertices[cur_edge.end]
        //);

        // Move to the end of this edge
        cur_vi = cur_edge.end;
        visited.insert(cur_vi);
        if cur_vi == start_vi {
            // We are done
            break;
        }
    }

    Ok(polygon)
}

#[test]
fn test_triangulate() {
    use ngeom::ops::*;
    use ngeom::re2::{AntiEven, AntiScalar, Vector};

    // TYPES SETUP

    struct Space2D;
    impl Space for Space2D {
        type Scalar = f32;
        type AntiScalar = AntiScalar<f32>;
        type Vector = Vector<f32>;
        type AxisOfRotation = Vector<f32>;
        type AntiEven = AntiEven<f32>;
    }

    #[derive(Clone)]
    struct VecVertex(Vec<Vector<f32>>);

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

    struct VecEdge(Vec<Edge<Space2D, usize>>);

    impl EdgeCollection for VecEdge {
        type Space = Space2D;
        type Index = usize;
        type VertexIndex = usize;

        fn iter(&self) -> impl Iterator<Item = &Edge<Space2D, usize>> {
            self.0.iter()
        }
        //fn iter_outgoing(&self, vi: usize) -> impl Iterator<Item = &Edge<Space2D, usize>> {
        //    self.0.iter().filter(move |e| e.start == vi)
        //}
    }

    struct VecFace(Vec<Face<usize>>);

    impl FaceCollection for VecFace {
        type Space = Space2D;
        type Index = usize;
        type EdgeIndex = usize;
        type VertexIndex = usize;

        fn iter(&self) -> impl Iterator<Item = &Face<usize>> {
            self.0.iter()
        }
        //fn iter_outgoing(&self, vi: usize) -> impl Iterator<Item = &Edge<Space2D, usize>> {
        //    self.0.iter().filter(move |e| e.start == vi)
        //}
    }
    // DATA

    let vertices = VecVertex(vec![
        Vector::point([0., 0.]),
        Vector::point([1., 0.]),
        Vector::point([0., 1.]),
    ]);

    let edges = VecEdge(vec![
        Edge::Line(Line { start: 0, end: 1 }),
        Edge::Arc(Arc {
            start: 1,
            axis: Vector::point([0., 0.]),
            end_angle: std::f32::consts::TAU / 8.,
            end: 2,
        }),
        Edge::Line(Line { start: 2, end: 3 }),
    ]);

    let boundary = FaceBoundary {
        edges: vec![(0, Dir::Fwd), (1, Dir::Fwd), (2, Dir::Fwd)],
    };

    let faces = VecFace(vec![Face {
        boundaries: vec![boundary],
    }]);

    println!("{:?}", triangulate(&vertices, &edges, &faces).unwrap());
}