use ngeom::algebraic_ops::*;
use ngeom::ops::*;
use ngeom::scalar::*;
use std::collections::BTreeSet;
use std::fmt::Debug;
use std::ops::{Index, Mul};

pub trait Space {
    type Scalar: Copy;
    type AntiScalar: Copy;
    type Vector: Copy + Transform<Self::AntiEven, Output = Self::Vector>;
    type AxisOfRotation: Copy
        + Mul<f32, Output = Self::AxisOfRotation>
        + Mul<Self::Scalar, Output = Self::AxisOfRotation>
        + AntiMul<Self::AntiScalar, Output = Self::AxisOfRotation>;
    type AntiEven: Copy + Rotor<Self::AxisOfRotation>;
}

pub trait VertexCollection: Index<Self::Index, Output = <Self::Space as Space>::Vector> {
    type Space: Space;
    type Index: Copy + Ord;
}

#[derive(Clone, Copy)]
pub struct Edge<SPACE: Space, VI: Copy + Ord> {
    pub start: VI,
    pub curve: Curve<SPACE>,
    pub end: VI,
}

pub trait EdgeCollection {
    type Space: Space;
    type Index: Copy;
    type VertexIndex: Copy + Ord;

    fn iter(&self) -> impl Iterator<Item = &Edge<Self::Space, Self::VertexIndex>>;
    fn iter_outgoing(
        &self,
        vi: Self::VertexIndex,
    ) -> impl Iterator<Item = &Edge<Self::Space, Self::VertexIndex>>;
}

#[derive(Clone, Copy, Debug)]
pub enum Curve<SPACE: Space> {
    Line,
    Arc(SPACE::AxisOfRotation),
    //CubicBezier--with optional offset
}

#[derive(Debug)]
pub enum PatchTopologyError {
    Empty,
    Open,
    Branching,
}

pub fn interpolate_patch<
    VC: VertexCollection,
    EC: EdgeCollection<Space = VC::Space, VertexIndex = VC::Index>,
>(
    vertices: &VC,
    edges: &EC,
) -> Result<Vec<<VC::Space as Space>::Vector>, PatchTopologyError>
where
    VC::Index: Debug,                    // TEMP
    <VC::Space as Space>::Vector: Debug, // TEMP
{
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

        // Interpolate this edge
        let start_pt = vertices[cur_edge.start];
        match cur_edge.curve {
            Curve::Line => {
                polygon.push(start_pt);
            }
            Curve::Arc(axis_of_rotation) => {
                for i in 0..10 {
                    // TODO dynamic spacing
                    let t = i as f32 / 10.;
                    polygon
                        .push(start_pt.transform(<VC::Space as Space>::AntiEven::rotor(
                            axis_of_rotation * t,
                        )));
                }
            }
        }
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
fn test_to_polygon() {
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

    struct VecVertex(Vec<Vector<f32>>);

    impl Index<usize> for VecVertex {
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
        fn iter_outgoing(&self, vi: usize) -> impl Iterator<Item = &Edge<Space2D, usize>> {
            self.0.iter().filter(move |e| e.start == vi)
        }
    }

    // DATA

    let vertices = VecVertex(vec![
        Vector::point([0., 0.]),
        Vector::point([1., 0.]),
        Vector::point([0., 1.]),
    ]);

    let edges = VecEdge(vec![
        Edge {
            start: 0,
            curve: Curve::Line,
            end: 1,
        },
        Edge {
            start: 1,
            curve: Curve::Arc(Vector::point([0., 0.]) * (std::f32::consts::TAU / 8.)),
            end: 2,
        },
        Edge {
            start: 2,
            curve: Curve::Line,
            end: 3,
        },
    ]);

    println!("{:?}", interpolate_patch(&vertices, &edges).unwrap());
}
