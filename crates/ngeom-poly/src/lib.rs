use enum_dispatch::enum_dispatch;
use ngeom::algebraic_ops::*;
use ngeom::ops::*;
use ngeom::scalar::*;
use std::cmp::Ordering;
use std::fmt;
use std::fmt::Debug;
use std::iter;
use std::ops::{Add, Index as StdIndex, Mul, Sub};

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

pub trait SpaceUv {
    type Scalar: Copy + Ring + Ord;
    type AntiScalar: Copy
        + Default
        + Ord
        + AntiMul<Self::AntiScalar, Output = Self::AntiScalar>
        + Mul<Self::Scalar, Output = Self::AntiScalar>;
    type Vector: Copy
        + Dot<Self::Vector, Output = Self::Scalar>
        + XHat
        + YHat
        + Join<Self::Vector, Output = Self::Bivector>
        + WeightNorm<Output = Self::AntiScalar>
        + BulkNorm<Output = Self::Scalar>
        + Unitized<Output = Self::Vector>
        + Sub<Self::Vector, Output = Self::Vector>;
    type Bivector: Copy
        + Join<Self::Vector, Output = Self::AntiScalar>
        + Meet<Self::Bivector, Output = Self::Vector>
        + WeightNorm<Output = Self::AntiScalar>
        + BulkNorm<Output = Self::Scalar>;
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

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct GraphNode {
    pub first_outgoing_edge: Option<usize>,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct GraphEdge {
    pub to: usize,
    pub next_outgoing_edge: Option<usize>,
}

#[derive(Clone, Default, PartialEq, Eq)]
pub struct Graph {
    pub nodes: Vec<GraphNode>,
    pub edges: Vec<GraphEdge>,
}

pub struct Successors<'graph> {
    graph: &'graph Graph,
    current_edge_index: Option<usize>,
}

impl<'graph> Iterator for Successors<'graph> {
    type Item = usize;

    fn next(&mut self) -> Option<usize> {
        match self.current_edge_index {
            None => None,
            Some(edge_index) => {
                let edge = &self.graph.edges[edge_index];
                self.current_edge_index = edge.next_outgoing_edge;
                Some(edge.to)
            }
        }
    }
}

pub struct Edges<'graph> {
    graph: &'graph Graph,
    current_node_index: usize,
    current_edge_index: Option<usize>,
}

impl<'graph> Iterator for Edges<'graph> {
    type Item = (usize, usize);

    fn next(&mut self) -> Option<(usize, usize)> {
        loop {
            if let Some(edge_index) = self.current_edge_index {
                let edge = &self.graph.edges[edge_index];
                // Advance iterator
                self.current_edge_index = edge.next_outgoing_edge;
                // Return this edge
                return Some((self.current_node_index, edge.to));
            } else {
                if self.current_node_index >= self.graph.nodes.len() - 1 {
                    return None; // This is the last node and it has no more edges, so iteration is finished
                } else {
                    // Move to the next node & loop
                    self.current_node_index += 1;
                    self.current_edge_index =
                        self.graph.nodes[self.current_node_index].first_outgoing_edge;
                }
            }
        }
    }
}

impl Graph {
    pub fn add_node(&mut self) -> usize {
        let index = self.nodes.len();
        self.nodes.push(Default::default());
        index
    }

    pub fn add_nodes(&mut self, count: usize) {
        self.nodes.reserve(count);
        self.nodes
            .extend(iter::repeat(Default::default()).take(count));
    }

    pub fn push_edge(&mut self, from: usize, to: usize) {
        let index = self.edges.len();
        let from_node = &mut self.nodes[from];
        self.edges.push(GraphEdge {
            to,
            next_outgoing_edge: from_node.first_outgoing_edge,
        });
        from_node.first_outgoing_edge = Some(index);
    }

    pub fn pop_edge(&mut self) -> Option<(usize, usize)> {
        for (from_index, from_node) in self.nodes.iter_mut().enumerate() {
            if let Some(edge_index) = from_node.first_outgoing_edge {
                let edge = &self.edges[edge_index];
                from_node.first_outgoing_edge = edge.next_outgoing_edge;
                return Some((from_index, edge.to));
            }
        }
        None
    }

    pub fn pop_outgoing_edge(&mut self, from: usize) -> Option<usize> {
        let from_node = &mut self.nodes[from];
        let edge = &self.edges[from_node.first_outgoing_edge?];
        from_node.first_outgoing_edge = edge.next_outgoing_edge;
        Some(edge.to)
    }

    pub fn remove_edge(&mut self, from: usize, to: usize) {
        let from_node = &mut self.nodes[from];

        let Some(edge_index) = from_node.first_outgoing_edge else {
            panic!("Requested node has no outgoing edges")
        };

        {
            let edge = &self.edges[edge_index];
            if edge.to == to {
                // Pop from node
                from_node.first_outgoing_edge = edge.next_outgoing_edge;
                return;
            }
        }

        let mut prev_edge_index = edge_index;

        loop {
            let Some(edge_index) = self.edges[prev_edge_index].next_outgoing_edge else {
                panic!("Requested edge not found")
            };

            let edge = &self.edges[edge_index];
            if edge.to == to {
                // Pop
                self.edges[prev_edge_index].next_outgoing_edge = edge.next_outgoing_edge;
                return;
            }
            prev_edge_index = edge_index;
        }
    }

    pub fn push_loop(&mut self, count: usize) {
        let start = self.nodes.len();
        let end = start + count;

        self.add_nodes(count);
        self.edges.reserve(count);

        match count {
            0 => {}
            1 => {
                self.push_edge(start, start);
            }
            _ => {
                for i in start..end - 1 {
                    self.push_edge(i, i + 1)
                }
                self.push_edge(end - 1, start);
            }
        }
    }

    pub fn new() -> Self {
        Self::default()
    }

    pub fn iter_successors(&self, from: usize) -> impl Iterator<Item = usize> + use<'_> {
        Successors {
            graph: self,
            current_edge_index: self.nodes[from].first_outgoing_edge,
        }
    }

    pub fn iter_edges(&self) -> impl Iterator<Item = (usize, usize)> + use<'_> {
        Edges {
            graph: self,
            current_node_index: 0,
            current_edge_index: self.nodes.get(0).and_then(|n| n.first_outgoing_edge),
        }
    }

    pub fn pop_min_outgoing_edge_by<F: Fn(&usize, &usize) -> Ordering>(
        &mut self,
        prev: usize,
        from: usize,
        cmp: F,
    ) -> Option<usize> {
        let to = self
            .iter_successors(from)
            .filter(|&n| n != prev)
            .min_by(cmp)?;
        self.remove_edge(from, to);

        Some(to)
    }
}

impl fmt::Debug for Graph {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Graph {{")?;
        let mut first_from = true;
        for from in 0..self.nodes.len() {
            if !first_from {
                write!(f, ", ")?;
            } else {
                first_from = false;
            }
            write!(f, "{:?} -> {{", from)?;
            let mut first_to = true;
            for to in self.iter_successors(from) {
                if !first_to {
                    write!(f, ", ")?;
                } else {
                    first_to = false;
                }
                write!(f, "{:?}", to)?;
            }
            write!(f, "}}")?;
        }
        write!(f, "}}")?;
        Ok(())
    }
}

pub struct Interpolation<POINT, UV> {
    pub points: Vec<POINT>,
    pub uv: Vec<UV>,
    pub graph: Graph,
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
    let mut graph = Graph::new();

    for face in faces.iter() {
        for boundary in face.boundaries.iter() {
            let start = points.len();
            for &(edge_index, dir) in boundary.edges.iter() {
                edges[edge_index].interpolate(vertices, dir, &mut points);
            }
            let end = points.len();
            graph.push_loop(end - start);
        }
    }

    // Now, re-interpret those boundary points into (U, V) coordinates.
    let uv: Vec<_> = points.iter().cloned().map(into_uv).collect();

    Interpolation { points, uv, graph }
}

fn u<VECTOR: Sized + Dot<VECTOR> + XHat>(pt: VECTOR) -> <VECTOR as Dot<VECTOR>>::Output {
    pt.dot(VECTOR::x_hat())
}
fn v<VECTOR: Sized + Dot<VECTOR> + YHat>(pt: VECTOR) -> <VECTOR as Dot<VECTOR>>::Output {
    pt.dot(VECTOR::y_hat())
}
fn cmp_uv<VECTOR: Sized + Copy + Dot<VECTOR, Output: Ord> + XHat + YHat>(
    uv1: VECTOR,
    uv2: VECTOR,
) -> Ordering {
    v(uv1).cmp(&v(uv2)).reverse().then(u(uv1).cmp(&u(uv2)))
}
fn is_on_left_chain<
    VECTOR: Copy + Join<VECTOR, Output: Join<VECTOR, Output: Default + PartialOrd>> + XHat + YHat,
>(
    pt: VECTOR,
    next_pt: VECTOR,
) -> bool {
    let anti_zero = Default::default();
    let going_down = pt.join(next_pt).join(VECTOR::x_hat());
    if going_down == anti_zero {
        // Edge is horizontal
        let going_right = pt.join(next_pt).join(VECTOR::y_hat());
        // Edge going right => left chain
        // Edge going left => right chain
        assert!(going_right != anti_zero, "Coincident points");
        going_right > anti_zero
    } else {
        // Edge going down => left chain
        // Edge going up => right chain
        going_down > anti_zero
    }
}

fn cmp_turn_angle<
    SCALAR: Copy + Default + Ord + Mul<SCALAR, Output = SCALAR>,
    VECTOR: Copy
        + Join<VECTOR, Output = BIVECTOR>
        + Sub<VECTOR, Output = VECTOR>
        + BulkNorm<Output = SCALAR>
        + Dot<VECTOR, Output = SCALAR>,
    BIVECTOR: Copy + BulkNorm<Output = SCALAR>,
    //ANTISCALAR: Copy
    //    + AntiMul<ANTISCALAR, Output = ANTISCALAR>
    //    + Default
    //    + Ord
    //    + Mul<SCALAR, Output = ANTISCALAR>,
>(
    pt1: VECTOR,
    pt2: VECTOR,
    pt3a: VECTOR,
    pt3b: VECTOR,
) -> Ordering {
    let v1 = pt2 - pt1;
    let v2a = pt3a - pt2;
    let v2b = pt3b - pt2;

    let a_w = v2a.bulk_norm();
    let b_w = v2b.bulk_norm();
    let a_sin = v1.join(v2a).bulk_norm();
    let b_sin = v1.join(v2b).bulk_norm();
    let a_cos = v1.dot(v2a);
    let b_cos = v1.dot(v2b);

    // For a valid comparison, we must normalize by the outgoing segment weight
    // which may be different between options A and B.
    // We can do this without division by moving the weight to the opposite side of the equation,
    // since both weights are positive.

    if a_sin >= Default::default() && b_sin >= Default::default() {
        // Both angles are not reflex,
        // so compare their cosines
        (b_w * a_cos).cmp(&(a_w * b_cos))
    } else if a_sin < Default::default() && b_sin < Default::default() {
        // Both angles are reflex
        // so compare their cosines
        (a_w * b_cos).cmp(&(b_w * a_cos))
    } else {
        // One angle is reflex and one is normal.
        // That means one sine is positive and one is non-negative.
        // This trivial comparison will return the appropriate result
        b_sin.cmp(&a_sin)
    }
}

pub fn partition_into_monotone_components<'a, SPACE: SpaceUv>(
    uv: &[SPACE::Vector],
    mut graph: Graph,
) -> Graph
where
    SPACE::Scalar: Debug,
    SPACE::Vector: Debug,
{
    // And finally, perform an argsort of the list
    // to get the sweep-line iteration order
    let mut sweep_line_order: Vec<_> = (0..uv.len()).collect();
    sweep_line_order.sort_by(|&i, &j| cmp_uv(uv[i], uv[j]));

    #[derive(Clone, Debug)]
    struct PolyEdge {
        from_node: usize,
        helper: usize,
        helper_is_merge: bool,
    }

    // The state of the sweep line algorithm: a sorted list of all active edges and their helpers.
    // (consider replacing this with a data structure allowing efficient insertion into the middle)
    let mut t = Vec::<PolyEdge>::new();

    fn search_t<
        VECTOR: Copy
            + Dot<VECTOR, Output: Copy + PartialOrd + Default>
            + Join<VECTOR, Output: Copy + Meet<<VECTOR as Join<VECTOR>>::Output, Output = VECTOR>>
            + XHat
            + WeightNorm<Output: PartialEq + Default>
            + Unitized<Output = VECTOR>,
    >(
        t: &Vec<PolyEdge>,
        uv: &[VECTOR],
        adjacent: &[(usize, usize)],
        pt: VECTOR,
    ) -> usize
    where
        <VECTOR as Dot<VECTOR>>::Output: Debug, // TEMPORARY
    {
        let pt_u = u(pt);
        let sweep_line = pt.join(VECTOR::x_hat());
        t.binary_search_by(|&PolyEdge { from_node, .. }| {
            let (_, to_node) = adjacent[from_node];
            let edge_line = uv[from_node].join(uv[to_node]);
            let intersection_pt = sweep_line.meet(edge_line);
            let intersection_u = u(if intersection_pt.weight_norm() == Default::default() {
                uv[to_node]
            } else {
                intersection_pt.unitized()
            });

            if intersection_u <= pt_u {
                Ordering::Less
            } else {
                Ordering::Greater
            }
        })
        .unwrap_err()
    }

    fn pop_t<
        VECTOR: Copy
            + Dot<VECTOR, Output: Copy + PartialOrd + Default>
            + Join<VECTOR, Output: Copy + Meet<<VECTOR as Join<VECTOR>>::Output, Output = VECTOR>>
            + XHat
            + WeightNorm<Output: PartialEq + Default>
            + Unitized<Output = VECTOR>,
    >(
        t: &mut Vec<PolyEdge>,
        uv: &[VECTOR],
        adjacent: &[(usize, usize)],
        pt: VECTOR,
    ) -> Option<PolyEdge>
    where
        <VECTOR as Dot<VECTOR>>::Output: Debug, // TEMPORARY
    {
        let ix = search_t(t, &uv, adjacent, pt);

        let ix = if ix == 0 { None } else { Some(ix - 1) }?;

        let edge = t.remove(ix);
        Some(edge)
    }

    fn push_t<
        VECTOR: Copy
            + Dot<VECTOR, Output: Copy + PartialOrd + Default>
            + Join<VECTOR, Output: Copy + Meet<<VECTOR as Join<VECTOR>>::Output, Output = VECTOR>>
            + XHat
            + WeightNorm<Output: PartialEq + Default>
            + Unitized<Output = VECTOR>,
    >(
        t: &mut Vec<PolyEdge>,
        uv: &[VECTOR],
        adjacent: &[(usize, usize)],
        pt: VECTOR,
        edge: PolyEdge,
    ) where
        <VECTOR as Dot<VECTOR>>::Output: Debug, // TEMPORARY
    {
        let i = search_t(t, &uv, adjacent, pt);
        t.insert(i, edge);
    }

    fn peek_t_mut<
        'a,
        VECTOR: Copy
            + Dot<VECTOR, Output: Copy + PartialOrd + Default>
            + Join<VECTOR, Output: Copy + Meet<<VECTOR as Join<VECTOR>>::Output, Output = VECTOR>>
            + XHat
            + WeightNorm<Output: PartialEq + Default>
            + Unitized<Output = VECTOR>,
    >(
        t: &'a mut Vec<PolyEdge>,
        uv: &[VECTOR],
        adjacent: &[(usize, usize)],
        pt: VECTOR,
    ) -> Option<&'a mut PolyEdge>
    where
        <VECTOR as Dot<VECTOR>>::Output: Debug, // TEMPORARY
    {
        let ix = search_t(t, &uv, adjacent, pt);

        let ix = if ix == 0 { None } else { Some(ix - 1) }?;

        Some(&mut t[ix])
    }

    /*
    let search_t =
        |t: &mut Vec<PolyEdge>, adjacent: &[(usize, usize)], pt: SPACE::Vector| -> usize {
            let pt_u = u(pt);
            let sweep_line = pt.join(SPACE::Vector::x_hat());
            t.binary_search_by(|&PolyEdge { from_node, .. }| {
                let (_, to_node) = adjacent[from_node];
                let edge_line = uv[from_node].join(uv[to_node]);
                let intersection_u = u(sweep_line.meet(edge_line).unitized());

                if intersection_u <= pt_u {
                    Ordering::Less
                } else {
                    Ordering::Greater
                }
            })
            .unwrap_err()
        };

    let pop_t = |t: &mut Vec<PolyEdge>,
                 adjacent: &[(usize, usize)],
                 pt: SPACE::Vector|
     -> Option<PolyEdge> {
        let ix = search_t(t, &uv, adjacent, pt);

        let ix = if ix == 0 { None } else { Some(ix - 1) }?;

        let edge = t.remove(ix);
        Some(edge)
    };

    let push_t =
        |t: &mut Vec<PolyEdge>, adjacent: &[(usize, usize)], pt: SPACE::Vector, edge: PolyEdge| {
            let i = search_t(t, &uv, adjacent, pt);
            t.insert(i, edge);
        };


    let peek_t_mut = |t: &'a mut Vec<PolyEdge>,
                      adjacent: &[(usize, usize)],
                      pt: SPACE::Vector|
     -> Option<&'a mut PolyEdge> {
        let ix = search_t(t, &uv, adjacent, pt);

        let ix = if ix == 0 { None } else { Some(ix - 1) }?;

        Some(&mut t[ix])
    };
    */

    //let pop_t_for_update = |t: &mut BTreeMap<SPACE::Scalar, PolyEdge>,
    //                        adjacent: &[(usize, usize)],
    //                        pt: SPACE::Vector|
    // -> Option<(SPACE::Scalar, PolyEdge)> {
    //    let edge = pop_t(t, pt)?;
    //    let from_node = edge.from_node;
    //    let (_, to_node) = adjacent[from_node];

    //    let sweep_line = pt.join(SPACE::Vector::x_hat());
    //    let edge_line = uv[from_node].join(uv[to_node]);

    //    // TODO this unitized() can be removed in favor of using a homogeneous magnitude
    //    let new_u = u(sweep_line.meet(edge_line).unitized());
    //    Some((new_u, edge))
    //};

    // Populate "next" and "prev" lists of nodes, from the graph
    let adjacent = {
        let empty = graph.nodes.len();
        let mut adjacent = vec![(empty, empty); graph.nodes.len()];
        for (
            i,
            GraphNode {
                first_outgoing_edge: mut cur_edge,
            },
        ) in graph.nodes.iter().enumerate()
        {
            while let Some(edge) = cur_edge {
                let from = i;
                let to = graph.edges[edge].to;
                assert!(adjacent[from].1 == empty, "Graph topology error--branching");
                assert!(adjacent[to].0 == empty, "Graph topology error--branching");
                adjacent[from].1 = to;
                adjacent[to].0 = from;

                cur_edge = graph.edges[edge].next_outgoing_edge;
            }
        }
        assert!(
            adjacent
                .iter()
                .all(|&(prev, next)| prev < empty && next < empty),
            "Graph topology error--dead ends"
        );
        adjacent
    };

    // Iterate over vertices in sweep-line order
    for i in sweep_line_order.into_iter() {
        // Figure out the type of vertex i by examining its neighbors
        let pt = uv[i];
        let (prev, next) = adjacent[i];
        let prev_pt = uv[prev];
        let next_pt = uv[next];
        println!("*** Point {:?} is {:?} ***", i, pt);
        println!("Prev is {:?}", prev_pt);
        println!("Next is {:?}", next_pt);
        println!("T is {:?}", t);

        let next_pt_below = match cmp_uv(next_pt, pt) {
            Ordering::Equal => panic!("Found coincident points"),
            Ordering::Greater => true,
            Ordering::Less => false,
        };
        let prev_pt_below = match cmp_uv(prev_pt, pt) {
            Ordering::Equal => panic!("Found coincident points"),
            Ordering::Greater => true,
            Ordering::Less => false,
        };

        // Calculate the direction we are turning at this vertex
        // to differentiate between start / split vertices and end / merge vertices
        let turn_direction = prev_pt.join(pt).join(next_pt);
        let anti_zero = SPACE::AntiScalar::default();

        if next_pt_below && prev_pt_below {
            if turn_direction >= anti_zero {
                // This is a start vertex
                println!("This is a start vertex");
                push_t(
                    &mut t,
                    &uv,
                    &adjacent,
                    pt,
                    PolyEdge {
                        from_node: i,
                        helper: i,
                        helper_is_merge: false,
                    },
                );
            } else {
                // This is a split vertex
                println!("This is a split vertex");
                {
                    let edge = peek_t_mut(&mut t, &uv, &adjacent, pt)
                        .expect("No active edge found to the left of point");
                    graph.push_edge(i, edge.helper);
                    graph.push_edge(edge.helper, i);
                    edge.helper = i;
                    edge.helper_is_merge = false;
                }
                push_t(
                    &mut t,
                    &uv,
                    &adjacent,
                    pt,
                    PolyEdge {
                        from_node: i,
                        helper: i,
                        helper_is_merge: false,
                    },
                );
            }
        } else if !next_pt_below && !prev_pt_below {
            if turn_direction > anti_zero {
                // This is a end vertex
                println!("This is a end vertex");
                let edge = pop_t(&mut t, &uv, &adjacent, pt)
                    .expect("No active edge found to the left of point");
                if edge.helper_is_merge {
                    graph.push_edge(i, edge.helper);
                    graph.push_edge(edge.helper, i);
                }
            } else {
                // This is a merge vertex
                println!("This is a merge vertex");
                let edge = pop_t(&mut t, &uv, &adjacent, pt)
                    .expect("No active edge found to the left of point");

                if edge.helper_is_merge {
                    graph.push_edge(i, edge.helper);
                    graph.push_edge(edge.helper, i);
                }

                let edge = peek_t_mut(&mut t, &uv, &adjacent, pt)
                    .expect("No active edge found to the left of point");
                if edge.helper_is_merge {
                    graph.push_edge(i, edge.helper);
                    graph.push_edge(edge.helper, i);
                }
                edge.helper = i;
                edge.helper_is_merge = true;
            }
        } else {
            // This is a normal vertex
            if is_on_left_chain(pt, next_pt) {
                // This vertex lies on the left side of the interval
                println!("This is a normal left vertex");
                let edge = pop_t(&mut t, &uv, &adjacent, pt)
                    .expect("No active edge found to the left of point");
                if edge.helper_is_merge {
                    // Walk the part of the polygon we are cutting off
                    graph.push_edge(i, edge.helper);
                    graph.push_edge(edge.helper, i);
                }
                push_t(
                    &mut t,
                    &uv,
                    &adjacent,
                    pt,
                    PolyEdge {
                        from_node: i,
                        helper: i,
                        helper_is_merge: false,
                    },
                );
            } else {
                // This vertex lies on the right side of the interval
                println!("This is a normal right vertex {:?}", pt);
                let edge = peek_t_mut(&mut t, &uv, &adjacent, pt)
                    .expect("No active edge found to the left of point");
                if edge.helper_is_merge {
                    graph.push_edge(i, edge.helper);
                    graph.push_edge(edge.helper, i);
                }
                edge.helper = i;
                edge.helper_is_merge = false;
            }
        }
    }

    graph
}

pub fn triangulate_monotone_components<SPACE: SpaceUv>(
    uv: &[SPACE::Vector],
    mut graph: Graph,
) -> Graph
where
    SPACE::Vector: Debug,     // TEMP
    SPACE::Scalar: Debug,     // TEMP
    SPACE::AntiScalar: Debug, // TEMP
{
    let anti_zero = SPACE::AntiScalar::default();

    // Initialize an "unvisited" set
    let mut unvisited_graph = graph.clone();
    let mut monotone_poly_edges = Vec::<(usize, usize)>::new();
    let mut s = Vec::<(usize, bool)>::new();

    while let Some((start_node, second_node)) = unvisited_graph.pop_edge() {
        // Walk a single monotone polygon from the graph
        monotone_poly_edges.clear(); // Prepare for a fresh polygon
        monotone_poly_edges.push((start_node, second_node));

        println!("START WALK {:?} to {:?}", start_node, second_node);

        let (mut prev_node, mut cur_node) = (start_node, second_node);

        loop {
            println!("walk from {:?} to {:?}", prev_node, cur_node);

            // Find the next node
            let uv_prev = uv[prev_node];
            let uv_cur = uv[cur_node];

            let next_node = unvisited_graph
                .pop_min_outgoing_edge_by(prev_node, cur_node, |&i, &j| {
                    cmp_turn_angle(uv_prev, uv_cur, uv[i], uv[j])
                })
                .expect("Bad manifold--no outgoing edge");
            monotone_poly_edges.push((cur_node, next_node));

            (prev_node, cur_node) = (cur_node, next_node);

            if cur_node == start_node {
                println!("DONE!");
                break;
            }
        }

        if monotone_poly_edges.len() > 3 {
            // Sort the current monotone polygon nodes in sweep-line order
            monotone_poly_edges.sort_by(|&(i, _), &(j, _)| cmp_uv(uv[i], uv[j]));

            // Initialize stack
            s.clear();
            s.push((monotone_poly_edges[0].0, false));
            s.push((
                monotone_poly_edges[1].0,
                is_on_left_chain(uv[monotone_poly_edges[1].0], uv[monotone_poly_edges[1].1]),
            ));
            println!(
                "Push {:?} {:?}",
                monotone_poly_edges[0], uv[monotone_poly_edges[0].0]
            );
            println!(
                "Push {:?} {:?} on the {}",
                monotone_poly_edges[1],
                uv[monotone_poly_edges[1].0],
                if is_on_left_chain(uv[monotone_poly_edges[1].0], uv[monotone_poly_edges[1].1]) {
                    "left"
                } else {
                    "right"
                }
            );
            for &(node, next) in &monotone_poly_edges[2..monotone_poly_edges.len() - 1] {
                let is_left = is_on_left_chain(uv[node], uv[next]);
                println!(
                    "Considering node {:?} {:?} on the {}",
                    node,
                    uv[node],
                    if is_left { "left" } else { "right" }
                );
                println!("Stack is {:?}", s);
                let (mut s_top, mut s_top_is_left) =
                    s.pop().expect("Stack empty during triangulation?");
                println!(
                    "Pop {:?} on the {}",
                    s_top,
                    if is_left { "left" } else { "right" }
                );
                if s_top_is_left != is_left {
                    // Different chains
                    println!("Different chain");
                    println!("Add diagonal {:?} to {:?}", node, s_top);
                    graph.push_edge(node, s_top);
                    graph.push_edge(s_top, node);
                    while s.len() > 1 {
                        let (s_top, _) = s.pop().unwrap();

                        // Add an edge to all nodes on the stack,
                        // except the last one
                        println!("Add diagonal {:?} to {:?}", node, s_top);
                        graph.push_edge(node, s_top);
                        graph.push_edge(s_top, node);
                    }
                    s.pop().expect("Stack empty during triangulation?"); // Discard last stack element
                } else {
                    // Same chain
                    println!("Same chain");
                    loop {
                        let Some(&(s_test_top, _)) = s.last() else {
                            break;
                        };

                        // See if adding an edge is possible
                        let proposed_tri = uv[s_test_top].join(uv[s_top]).join(uv[node]);
                        if match is_left {
                            true => proposed_tri < anti_zero,
                            false => proposed_tri > anti_zero,
                        } {
                            println!(
                                "Diag {:?} to {:?} doesn't work-- tri={:?}-{:?}-{:?}={:?}",
                                node, s_test_top, s_test_top, s_top, node, proposed_tri
                            );
                            // This diagonal doesn't work
                            break;
                        }

                        (s_top, s_top_is_left) = s.pop().unwrap();

                        println!("Add diagonal {:?} to {:?}", node, s_top);
                        graph.push_edge(node, s_top);
                        graph.push_edge(s_top, node);
                    }
                }
                s.push((s_top, s_top_is_left));
                s.push((node, is_left));
            }

            // Last node: Add edges to all remaining nodes in the stack,
            // except the last one
            let &(node, _) = monotone_poly_edges.last().unwrap();
            println!("Handling last node: {:?} {:?}", node, uv[node]);
            println!("Stack is {:?}", s);
            for &(s_entry, _) in &s[1..s.len() - 1] {
                println!("Add diagonal {:?} to {:?}", node, s_entry);
                graph.push_edge(node, s_entry);
                graph.push_edge(s_entry, node);
            }
        }
    }

    graph
}

pub fn graph_to_triangles<SPACE: SpaceUv>(
    uv: &[SPACE::Vector],
    mut graph: Graph,
) -> Vec<[usize; 3]> {
    // Initialize an "unvisited" set
    let mut triangles = Vec::<[usize; 3]>::new();

    while let Some((node1, node2)) = graph.pop_edge() {
        // Walk a single monotone polygon from the graph

        println!("*TRI*");
        println!("Node 1: {:?}", node1);
        println!("Node 2: {:?}", node2);
        assert!(node2 != node1, "Graph has a 1-cycle");

        let node3 = {
            let uv_prev = uv[node1];
            let uv_cur = uv[node2];

            let node3 = graph
                .pop_min_outgoing_edge_by(node1, node2, |&i, &j| {
                    let result = cmp_turn_angle(uv_prev, uv_cur, uv[i], uv[j]);
                    println!(
                        "CMP: {:?}-{:?}-{:?} {:?} {:?}-{:?}-{:?}",
                        node1, node2, i, result, node1, node2, j
                    );
                    result
                })
                .expect("Bad manifold--no outgoing edge");
            assert!(node3 != node2, "Graph has a 1-cycle");
            assert!(node3 != node1, "Graph has a 2-cycle");
            println!("Node 3: {:?}", node3);
            node3
        };

        {
            let uv_prev = uv[node2];
            let uv_cur = uv[node3];

            let node4 = graph
                .pop_min_outgoing_edge_by(node2, node3, |&i, &j| {
                    let result = cmp_turn_angle(uv_prev, uv_cur, uv[i], uv[j]);
                    println!(
                        "CMP: {:?}-{:?}-{:?} {:?} {:?}-{:?}-{:?}",
                        node2, node3, i, result, node1, node2, j
                    );
                    result
                })
                .expect("Bad manifold--no outgoing edge");
            println!("Node 4: {:?}", node4);
            assert!(
                node4 == node1,
                "Graph has a >3 cycle (not fully triangulated)"
            );
        }

        println!("OK");

        triangles.push([node1, node2, node3]);
    }
    triangles
}

#[cfg(test)]
mod tests {
    use derive_more::*;
    use ngeom::ops::*;
    use ngeom::re2::{AntiEven, AntiScalar, Bivector, Vector};
    use ngeom::scalar::*;
    use std::cmp::Ordering;
    use std::ops::Index as StdIndex;

    use super::*;

    // TYPES SETUP

    #[derive(
        Clone,
        Copy,
        Default,
        Display,
        derive_more::Debug,
        Neg,
        Add,
        Sub,
        Mul,
        Div,
        AddAssign,
        SubAssign,
        MulAssign,
        DivAssign,
        PartialOrd,
        From,
    )]
    #[mul(forward)]
    #[div(forward)]
    #[debug("{:?}", self.0)]
    struct R32(f32);

    impl std::ops::Mul<f32> for R32 {
        type Output = R32;
        fn mul(self, r: f32) -> R32 {
            R32(self.0 * r)
        }
    }

    impl std::ops::Div<f32> for R32 {
        type Output = R32;
        fn div(self, r: f32) -> R32 {
            R32(self.0 * r)
        }
    }

    impl Abs for R32 {
        type Output = R32;
        fn abs(self) -> R32 {
            R32(self.0.abs())
        }
    }

    impl Ring for R32 {
        fn one() -> R32 {
            R32(1.)
        }
    }

    impl Rational for R32 {
        fn one_half() -> R32 {
            R32(1. / 2.)
        }
        fn one_third() -> R32 {
            R32(1. / 3.)
        }
        fn one_fifth() -> R32 {
            R32(1. / 5.)
        }
    }

    impl Sqrt for R32 {
        type Output = R32;
        fn sqrt(self) -> R32 {
            R32(self.0.sqrt())
        }
    }

    impl Trig for R32 {
        type Output = R32;

        fn cos(self) -> R32 {
            R32(self.0.cos())
        }
        fn sin(self) -> R32 {
            R32(self.0.sin())
        }
        fn sinc(self) -> R32 {
            let self_adj = self.0.abs() + f32::EPSILON;
            R32(self_adj.sin() / self_adj)
        }
    }

    impl Recip for R32 {
        type Output = R32;

        // This is not NaN-free, even for internal usage!
        fn recip(self) -> R32 {
            R32(self.0.recip())
        }
    }

    impl PartialEq for R32 {
        fn eq(&self, other: &R32) -> bool {
            if self.0.is_finite() && other.0.is_finite() {
                self.0 == other.0
            } else {
                panic!("Uncomparable float values");
            }
        }
    }
    impl Eq for R32 {}

    impl Ord for R32 {
        fn cmp(&self, r: &R32) -> Ordering {
            self.partial_cmp(r).expect("Uncomparable float values")
        }
    }

    struct Space2D;
    impl Space for Space2D {
        type Scalar = R32;
        type AntiScalar = AntiScalar<R32>;
        type Vector = Vector<R32>;
        type AxisOfRotation = Vector<R32>;
        type AntiEven = AntiEven<R32>;
    }
    impl SpaceUv for Space2D {
        type Scalar = R32;
        type AntiScalar = AntiScalar<R32>;
        type Vector = Vector<R32>;
        type Bivector = Bivector<R32>;
    }

    #[derive(Clone)]
    struct VecVertex(Vec<Vector<R32>>);

    impl StdIndex<usize> for VecVertex {
        type Output = Vector<R32>;

        fn index(&self, idx: usize) -> &Vector<R32> {
            self.0.index(idx)
        }
    }

    impl VertexCollection for VecVertex {
        type Space = Space2D;
        type Index = usize;
    }

    struct VecEdge(Vec<Edge<Space2D, usize>>);

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

    #[test]
    fn test_make_monotone_already_monotone() {
        {
            // Triangle
            let uv = vec![
                Vector::point([R32(0.), R32(0.)]),
                Vector::point([R32(1.), R32(0.)]),
                Vector::point([R32(0.), R32(1.)]),
            ];

            let mut graph = Graph::new();
            graph.push_loop(3);

            let new_graph = partition_into_monotone_components::<Space2D>(&uv, graph.clone());
            assert_eq!(graph, new_graph);
        }

        {
            // Square
            let uv = vec![
                Vector::point([R32(0.), R32(0.)]),
                Vector::point([R32(1.), R32(0.)]),
                Vector::point([R32(1.), R32(1.)]),
                Vector::point([R32(0.), R32(1.)]),
            ];

            let mut graph = Graph::new();
            graph.push_loop(4);

            let new_graph = partition_into_monotone_components::<Space2D>(&uv, graph.clone());
            assert_eq!(graph, new_graph);
        }

        {
            // Square with redundant edge points
            let uv = vec![
                Vector::point([R32(0.), R32(0.)]),
                Vector::point([R32(0.5), R32(0.)]),
                Vector::point([R32(1.), R32(0.)]),
                Vector::point([R32(1.), R32(0.5)]),
                Vector::point([R32(1.), R32(1.)]),
                Vector::point([R32(0.5), R32(1.)]),
                Vector::point([R32(0.), R32(1.)]),
                Vector::point([R32(0.), R32(0.5)]),
            ];

            let mut graph = Graph::new();
            graph.push_loop(8);

            let new_graph = partition_into_monotone_components::<Space2D>(&uv, graph.clone());
            assert_eq!(graph, new_graph);
        }
    }

    #[test]
    fn test_make_monotone_comb() {
        {
            // Three-pointed comb like this: /\/\/\
            let uv = vec![
                Vector::point([R32(6.), R32(0.)]),
                Vector::point([R32(5.), R32(1.)]),
                Vector::point([R32(4.), R32(0.5)]),
                Vector::point([R32(3.), R32(1.)]),
                Vector::point([R32(2.), R32(0.5)]),
                Vector::point([R32(1.), R32(1.)]),
                Vector::point([R32(0.), R32(0.)]),
            ];

            let mut graph = Graph::new();
            graph.push_loop(7);

            let graph = partition_into_monotone_components::<Space2D>(&uv, graph);
            // Two diagonals should have been added for a total of 4 more graph edges
            assert_eq!(graph.edges.len(), 11);
        }

        {
            // Upside-down three-pointed comb like this: \/\/\/
            let uv = vec![
                Vector::point([R32(0.), R32(0.)]),
                Vector::point([R32(1.), R32(-1.)]),
                Vector::point([R32(2.), R32(-0.5)]),
                Vector::point([R32(3.), R32(-1.)]),
                Vector::point([R32(4.), R32(-0.5)]),
                Vector::point([R32(5.), R32(-1.)]),
                Vector::point([R32(6.), R32(0.)]),
            ];

            let mut graph = Graph::new();
            graph.push_loop(7);

            let graph = partition_into_monotone_components::<Space2D>(&uv, graph);
            // Two diagonals should have been added for a total of 4 more graph edges
            assert_eq!(graph.edges.len(), 11);
        }
    }

    #[test]
    fn test_make_monotone_square_hole() {
        {
            // Square
            let uv = vec![
                // Outside
                Vector::point([R32(0.), R32(0.)]),
                Vector::point([R32(1.), R32(0.)]),
                Vector::point([R32(1.), R32(1.)]),
                Vector::point([R32(0.), R32(1.)]),
                // Hole
                Vector::point([R32(0.25), R32(0.75)]),
                Vector::point([R32(0.75), R32(0.75)]),
                Vector::point([R32(0.75), R32(0.25)]),
                Vector::point([R32(0.25), R32(0.25)]),
            ];

            let mut graph = Graph::new();
            graph.push_loop(4);
            graph.push_loop(4);

            let graph = partition_into_monotone_components::<Space2D>(&uv, graph);
            // Two diagonals should have been added for a total of 4 more graph edges
            assert_eq!(graph.edges.len(), 12);
        }
    }

    #[test]
    fn test_make_monotone_multiple_shapes() {
        // Square
        let uv = vec![
            // Outside
            Vector::point([R32(0.), R32(0.)]),
            Vector::point([R32(1.), R32(0.)]),
            Vector::point([R32(1.), R32(1.)]),
            Vector::point([R32(0.), R32(1.)]),
            // Hole
            Vector::point([R32(0.25), R32(0.75)]),
            Vector::point([R32(0.75), R32(0.75)]),
            Vector::point([R32(0.75), R32(0.25)]),
            Vector::point([R32(0.25), R32(0.25)]),
            // Second square outside
            Vector::point([R32(2.), R32(0.)]),
            Vector::point([R32(3.), R32(0.)]),
            Vector::point([R32(3.), R32(1.)]),
            Vector::point([R32(2.), R32(1.)]),
            // Second square hole
            Vector::point([R32(2.25), R32(0.75)]),
            Vector::point([R32(2.75), R32(0.75)]),
            Vector::point([R32(2.75), R32(0.25)]),
            Vector::point([R32(2.25), R32(0.25)]),
        ];

        let mut graph = Graph::new();
        graph.push_loop(4);
        graph.push_loop(4);
        graph.push_loop(4);
        graph.push_loop(4);

        let graph = partition_into_monotone_components::<Space2D>(&uv, graph);
        // Four diagonals should have been added for a total of 8 more graph edges
        assert_eq!(graph.edges.len(), 24);

        // TODO make sure no diagonals go from the first shape (vertices 0-8) to the second shape (vertices 8-16)
    }

    #[test]
    fn test_triangulate_simple_shapes() {
        {
            // Triangle
            let uv = vec![
                Vector::point([R32(0.), R32(0.)]),
                Vector::point([R32(1.), R32(0.)]),
                Vector::point([R32(0.), R32(1.)]),
            ];

            let mut graph = Graph::new();
            graph.push_loop(3);

            let graph = triangulate_monotone_components::<Space2D>(&uv, graph);
            let tri = graph_to_triangles::<Space2D>(&uv, graph);
            assert_eq!(tri.len(), 1);
        }

        {
            // Square
            let uv = vec![
                Vector::point([R32(0.), R32(0.)]),
                Vector::point([R32(1.), R32(0.)]),
                Vector::point([R32(1.), R32(1.)]),
                Vector::point([R32(0.), R32(1.)]),
            ];

            let mut graph = Graph::new();
            graph.push_loop(4);

            let graph = triangulate_monotone_components::<Space2D>(&uv, graph);
            println!("TRIANGULATE: {:?}", graph);
            println!("{:?}", uv);
            let tri = graph_to_triangles::<Space2D>(&uv, graph);
            assert_eq!(tri.len(), 2);
        }

        {
            // Square with redundant edge points
            let uv = vec![
                Vector::point([R32(0.), R32(0.)]),
                Vector::point([R32(0.5), R32(0.)]),
                Vector::point([R32(1.), R32(0.)]),
                Vector::point([R32(1.), R32(0.5)]),
                Vector::point([R32(1.), R32(1.)]),
                Vector::point([R32(0.5), R32(1.)]),
                Vector::point([R32(0.), R32(1.)]),
                Vector::point([R32(0.), R32(0.5)]),
            ];

            let mut graph = Graph::new();
            graph.push_loop(8);

            let graph = triangulate_monotone_components::<Space2D>(&uv, graph);
            println!("TRIANGULATE: {:?}", graph);
            println!("{:?}", uv);
            let tri = graph_to_triangles::<Space2D>(&uv, graph);
            assert_eq!(tri.len(), 6);
        }
    }

    #[test]
    fn test_triangulate_square_hole() {
        {
            // Square
            let uv = vec![
                // Outside
                Vector::point([R32(0.), R32(0.)]),
                Vector::point([R32(1.), R32(0.)]),
                Vector::point([R32(1.), R32(1.)]),
                Vector::point([R32(0.), R32(1.)]),
                // Hole
                Vector::point([R32(0.25), R32(0.75)]),
                Vector::point([R32(0.75), R32(0.75)]),
                Vector::point([R32(0.75), R32(0.25)]),
                Vector::point([R32(0.25), R32(0.25)]),
            ];

            let mut graph = Graph::new();
            graph.push_loop(4);
            graph.push_loop(4);

            let graph = partition_into_monotone_components::<Space2D>(&uv, graph);
            println!("Graph after partitioning into monotone is {:?}", graph);
            let graph = triangulate_monotone_components::<Space2D>(&uv, graph);
            println!("Graph after triangulation is {:?}", graph);
            let _tri = graph_to_triangles::<Space2D>(&uv, graph);
        }
    }

    //#[test]
    //fn test_triangulate() {
    //    // DATA
    //    let vertices = VecVertex(vec![
    //        Vector::point([R32(0.), R32(0.)]),
    //        Vector::point([R32(1.), R32(0.)]),
    //        Vector::point([R32(0.), R32(1.)]),
    //    ]);

    //    let edges = VecEdge(vec![
    //        Edge::Line(Line { start: 0, end: 1 }),
    //        Edge::Arc(Arc {
    //            start: 1,
    //            axis: Vector::point([R32(0.), R32(0.)]),
    //            end_angle: R32(std::f32::consts::TAU / 8.),
    //            end: 2,
    //        }),
    //        //Edge::Line(Line { start: 1, end: 2 }),
    //        Edge::Line(Line { start: 2, end: 3 }),
    //    ]);

    //    let boundary = FaceBoundary {
    //        edges: vec![(0, Dir::Fwd), (1, Dir::Fwd), (2, Dir::Fwd)],
    //    };

    //    let faces = VecFace(vec![Face {
    //        boundaries: vec![boundary],
    //    }]);

    //    let Interpolation {
    //        points: _,
    //        uv,
    //        graph,
    //    } = interpolate(&vertices, &edges, &faces, |pt| pt);

    //    let graph = partition_into_monotone_components::<Space2D>(&uv, graph);
    //    let graph = triangulate_monotone_components::<Space2D>(&uv, graph);
    //    let triangles = graph_to_triangles::<Space2D>(&uv, graph);
    //    println!("{:?}", triangles);
    //    //triangulate::<Space2D, _, _, _, _>(&vertices, &edges, &faces, |pt| pt);
    //    //println!("{:?}", triangulate(&vertices, &edges, &faces).unwrap());
    //}
}
