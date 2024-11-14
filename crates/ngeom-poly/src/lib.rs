use enum_dispatch::enum_dispatch;
use ngeom::algebraic_ops::*;
use ngeom::ops::*;
use ngeom::scalar::*;
use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Debug;
use std::iter;
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

pub trait SpaceUv {
    type Scalar: Copy + Ring + Ord;
    type AntiScalar: Copy + Default + PartialOrd;
    type Vector: Copy
        + Dot<Self::Vector, Output = Self::Scalar>
        + XHat
        + YHat
        + Join<Self::Vector, Output = Self::Bivector>
        + Unitized<Output = Self::Vector>;
    type Bivector: Copy
        + Join<Self::Vector, Output = Self::AntiScalar>
        + Meet<Self::Bivector, Output = Self::Vector>;
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

#[derive(Clone, Debug, Default)]
pub struct GraphNode {
    pub first_outgoing_edge: Option<usize>,
}

#[derive(Clone, Debug, Default)]
pub struct GraphEdge {
    pub to: usize,
    pub next_outgoing_edge: Option<usize>,
}

#[derive(Clone, Debug, Default)]
pub struct Graph {
    pub nodes: Vec<GraphNode>,
    pub edges: Vec<GraphEdge>,
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

    pub fn pop_edge(&mut self, from: usize) -> Option<usize> {
        let from_node = &mut self.nodes[from];
        let removed_edge = from_node.first_outgoing_edge;
        if let Some(removed_edge) = removed_edge {
            from_node.first_outgoing_edge = self.edges[removed_edge].next_outgoing_edge;
        }
        removed_edge
    }

    pub fn new() -> Self {
        Self::default()
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
            graph.add_nodes(end - start);
            if end - start == 0 {
                // Interpolation contained no points--mistake?
            } else if end - start == 1 {
                // Interpolation contained one point--include this point in the triangulation
                graph.push_edge(start, start);
                panic!("Cannot handle 1-point polygon (TODO)");
            } else {
                graph.edges.reserve(end - start);

                for i in start..end - 1 {
                    graph.push_edge(i, i + 1)
                }
                graph.push_edge(end - 1, start);

                //if end - start == 2 {
                //    panic!("Cannot handle 2-point polygon (TODO)");
                //}
            }
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
    VECTOR: Join<VECTOR, Output: Join<VECTOR, Output: Default + PartialOrd>> + XHat,
>(
    pt: VECTOR,
    next_pt: VECTOR,
) -> bool {
    let anti_zero = Default::default();
    let turn_direction = pt.join(next_pt).join(VECTOR::x_hat());
    turn_direction > anti_zero
}

pub fn partition_into_monotone_components<SPACE: SpaceUv>(
    uv: &[SPACE::Vector],
    mut graph: Graph,
) -> Graph
where
    SPACE::Vector: Debug,
    SPACE::Scalar: Debug,
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

    // The state of the sweep line algorithm: a binary tree containing edges and their helpers.
    // The key is the intersection of the edge and the sweep line
    let mut t = BTreeMap::<SPACE::Scalar, PolyEdge>::new();

    let pop_t =
        |t: &mut BTreeMap<SPACE::Scalar, PolyEdge>, pt: SPACE::Vector| -> Option<PolyEdge> {
            let (&old_u, _) = t.range(..=u(pt)).next_back()?;
            let edge = t.remove(&old_u).unwrap();
            Some(edge)
        };

    let pop_t_for_update = |t: &mut BTreeMap<SPACE::Scalar, PolyEdge>,
                            adjacent: &[(usize, usize)],
                            pt: SPACE::Vector|
     -> Option<(SPACE::Scalar, PolyEdge)> {
        let edge = pop_t(t, pt)?;
        let from_node = edge.from_node;
        let (_, to_node) = adjacent[from_node];

        let sweep_line = pt.join(SPACE::Vector::x_hat());
        let edge_line = uv[from_node].join(uv[to_node]);

        let new_u = u(sweep_line.meet(edge_line).unitized());
        Some((new_u, edge))
    };

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
        //println!("Point is {:?}", pt);
        //println!("Prev is {:?}", prev_pt);
        //println!("Next is {:?}", next_pt);
        //println!("T is {:?}", t);

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
                //println!("This is a start vertex");
                if let Some((k, v)) = pop_t_for_update(&mut t, &adjacent, pt) {
                    t.insert(k, v);
                }
                t.insert(
                    u(pt),
                    PolyEdge {
                        from_node: i,
                        helper: i,
                        helper_is_merge: false,
                    },
                );
            } else {
                // This is a split vertex
                //println!("This is a split vertex");
                let edge = pop_t(&mut t, pt).expect("No active edge found to the left of point");
                // TODO ADD DIAGONAL
                //println!("Draw diagonal {:?} to {:?}", i, edge.helper);
                t.insert(
                    u(pt),
                    PolyEdge {
                        from_node: edge.from_node,
                        helper: i,
                        helper_is_merge: false,
                    },
                );
            }
        } else if !next_pt_below && !prev_pt_below {
            if turn_direction > anti_zero {
                // This is a end vertex
                //println!("This is a end vertex");
                let edge = pop_t(&mut t, pt).expect("No active edge found to the left of point");
                if edge.helper_is_merge {
                    graph.push_edge(i, edge.helper);
                    graph.push_edge(edge.helper, i);
                }
            } else {
                // This is a merge vertex
                //println!("This is a merge vertex");
                let edge = pop_t(&mut t, pt).expect("No active edge found to the left of point");

                if edge.helper_is_merge {
                    graph.push_edge(i, edge.helper);
                    graph.push_edge(edge.helper, i);
                }

                let (k, mut edge) = pop_t_for_update(&mut t, &adjacent, pt)
                    .expect("No active edge found to the left of point");
                if edge.helper_is_merge {
                    graph.push_edge(i, edge.helper);
                    graph.push_edge(edge.helper, i);
                }
                edge.helper = i;
                edge.helper_is_merge = true;
                t.insert(k, edge);
            }
        } else {
            // This is a normal vertex
            if is_on_left_chain(pt, next_pt) {
                // This vertex lies on the left side of the interval
                //println!("This is a normal left vertex");
                let edge = pop_t(&mut t, pt).expect("No active edge found to the left of point");
                if edge.helper_is_merge {
                    // Walk the part of the polygon we are cutting off
                    graph.push_edge(i, edge.helper);
                    graph.push_edge(edge.helper, i);
                }
                t.insert(
                    u(pt),
                    PolyEdge {
                        from_node: i,
                        helper: i,
                        helper_is_merge: false,
                    },
                );
            } else {
                // This vertex lies on the right side of the interval
                //println!("This is a normal right vertex {:?}", pt);
                let (k, mut edge) = pop_t_for_update(&mut t, &adjacent, pt)
                    .expect("No active edge found to the left of point");
                if edge.helper_is_merge {
                    graph.push_edge(i, edge.helper);
                    graph.push_edge(edge.helper, i);
                }
                edge.helper = i;
                edge.helper_is_merge = false;
                t.insert(k, edge);
            }
        }
    }

    graph
}

pub fn triangulate_monotone_components<SPACE: SpaceUv>(
    uv: &[SPACE::Vector],
    mut graph: Graph,
) -> Graph {
    let anti_zero = SPACE::AntiScalar::default();

    // Initialize an "unvisited" set
    let mut unvisited_nodes = BTreeSet::<usize>::from_iter(0..graph.nodes.len());
    let mut unvisited_graph = graph.clone();
    let mut monotone_poly_nodes = Vec::<usize>::new();
    let mut s = Vec::<(usize, bool)>::new();

    while let Some(start_node) = unvisited_nodes.pop_first() {
        // Walk a single monotone polygon from the graph
        monotone_poly_nodes.clear(); // Prepare for a fresh polygon
        let mut cur_node = start_node;

        loop {
            monotone_poly_nodes.push(cur_node);

            // Find the next node
            let next_edge = unvisited_graph
                .pop_edge(cur_node)
                .expect("Bad manifold--no outgoing edge");
            cur_node = graph.edges[next_edge].to;

            if cur_node == start_node {
                break;
            }
            unvisited_nodes.remove(&cur_node);
        }

        if monotone_poly_nodes.len() > 3 {
            // Store "next" pointers for the current monotone polygon
            let next: Vec<_> = monotone_poly_nodes[1..]
                .iter()
                .cloned()
                .chain(Some(monotone_poly_nodes[0]).into_iter())
                .collect();

            // Sort the current monotone polygon nodes in sweep-line order
            monotone_poly_nodes.sort_by(|&i, &j| cmp_uv(uv[i], uv[j]));

            // Initialize stack
            s.clear();
            s.push((monotone_poly_nodes[0], false));
            s.push((
                monotone_poly_nodes[1],
                is_on_left_chain(uv[monotone_poly_nodes[1]], uv[next[monotone_poly_nodes[1]]]),
            ));
            for &node in &monotone_poly_nodes[2..monotone_poly_nodes.len() - 1] {
                let is_left = is_on_left_chain(uv[node], uv[next[node]]);
                let (mut s_top, mut s_top_is_left) =
                    s.pop().expect("Stack empty during triangulaton?");
                if s_top_is_left != is_left {
                    // Different chains
                    while !s.is_empty() {
                        // Add an edge to all nodes on the stack,
                        // except the last one
                        graph.push_edge(node, s_top);
                        graph.push_edge(s_top, node);
                        (s_top, s_top_is_left) = s.pop().unwrap();
                    }
                } else {
                    // Same chain
                    loop {
                        let Some((s_penultimate, s_penultimate_is_left)) = s.pop() else {
                            break;
                        };
                        // Add an edge if possible
                        let proposed_tri = uv[node].join(uv[s_top]).join(uv[s_penultimate]);
                        if is_left && proposed_tri > anti_zero
                            || !is_left && proposed_tri < anti_zero
                        {
                            break;
                        }

                        (s_top, s_top_is_left) = (s_penultimate, s_penultimate_is_left);

                        graph.push_edge(node, s_top);
                        graph.push_edge(s_top, node);
                    }
                }
                s.push((s_top, s_top_is_left));
                s.push((node, is_left));
            }
        }
    }

    graph
}

pub fn graph_to_triangles(mut graph: Graph) -> Vec<[usize; 3]> {
    // Initialize an "unvisited" set
    let mut unvisited_nodes = BTreeSet::<usize>::from_iter(0..graph.nodes.len());
    let mut triangles = Vec::<[usize; 3]>::new();

    while let Some(node1) = unvisited_nodes.pop_first() {
        // Walk a single monotone polygon from the graph

        let next_edge = graph
            .pop_edge(node1)
            .expect("Bad manifold--no outgoing edge");
        let node2 = graph.edges[next_edge].to;
        assert!(node2 != node1, "Graph has a 1-cycle");

        let next_edge = graph
            .pop_edge(node2)
            .expect("Bad manifold--no outgoing edge");
        let node3 = graph.edges[next_edge].to;
        assert!(node3 != node1, "Graph has a 2-cycle");

        let next_edge = graph
            .pop_edge(node3)
            .expect("Bad manifold--no outgoing edge");
        assert!(
            node1 == graph.edges[next_edge].to,
            "Graph has a >3 cycle (not fully triangulated)"
        );

        unvisited_nodes.remove(&node1);
        unvisited_nodes.remove(&node2);
        unvisited_nodes.remove(&node3);

        triangles.push([node1, node2, node3]);
    }
    triangles
}

#[test]
fn test_triangulate() {
    use derive_more::*;
    use ngeom::ops::*;
    use ngeom::re2::{AntiEven, AntiScalar, Bivector, Vector};

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
    // DATA

    let vertices = VecVertex(vec![
        Vector::point([R32(0.), R32(0.)]),
        Vector::point([R32(1.), R32(0.)]),
        Vector::point([R32(0.), R32(1.)]),
    ]);

    let edges = VecEdge(vec![
        Edge::Line(Line { start: 0, end: 1 }),
        //Edge::Arc(Arc {
        //    start: 1,
        //    axis: Vector::point([R32(0.), R32(0.)]),
        //    end_angle: R32(std::f32::consts::TAU / 8.),
        //    end: 2,
        //}),
        Edge::Line(Line { start: 1, end: 2 }),
        Edge::Line(Line { start: 2, end: 3 }),
    ]);

    let boundary = FaceBoundary {
        edges: vec![(0, Dir::Fwd), (1, Dir::Fwd), (2, Dir::Fwd)],
    };

    let faces = VecFace(vec![Face {
        boundaries: vec![boundary],
    }]);

    let Interpolation {
        points: _,
        uv,
        graph,
    } = interpolate(&vertices, &edges, &faces, |pt| pt);

    let graph = partition_into_monotone_components::<Space2D>(&uv, graph);
    let graph = triangulate_monotone_components::<Space2D>(&uv, graph);
    let triangles = graph_to_triangles(graph);
    println!("{:?}", triangles);
    //triangulate::<Space2D, _, _, _, _>(&vertices, &edges, &faces, |pt| pt);
    //println!("{:?}", triangulate(&vertices, &edges, &faces).unwrap());
}
