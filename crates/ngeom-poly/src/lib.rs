use enum_dispatch::enum_dispatch;
use ngeom::algebraic_ops::*;
use ngeom::ops::*;
use ngeom::scalar::*;
use std::cmp::Ordering;
use std::collections::BTreeSet;
use std::fmt;
use std::fmt::Debug;
use std::iter::Iterator;
use std::mem::replace;
use std::ops::{Add, Bound, Index as StdIndex, Mul, RangeBounds, Sub};

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

#[derive(Clone, Default, PartialEq, Eq)]
pub struct EdgeSet(BTreeSet<(usize, usize)>);

impl EdgeSet {
    pub fn new() -> Self {
        EdgeSet(BTreeSet::new())
    }
    pub fn insert(&mut self, from: usize, to: usize) -> bool {
        self.0.insert((from, to))
    }
    pub fn remove(&mut self, from: usize, to: usize) -> bool {
        self.0.remove(&(from, to))
    }
    pub fn pop(&mut self) -> Option<(usize, usize)> {
        self.0.pop_last()
    }
    pub fn insert_loop(&mut self, range: impl RangeBounds<usize>) {
        let start = match range.start_bound() {
            Bound::Included(&value) => value,
            Bound::Excluded(&value) => value + 1,
            Bound::Unbounded => panic!("Cannot add unbounded loop"),
        };
        let end = match range.end_bound() {
            Bound::Included(&value) => value + 1,
            Bound::Excluded(&value) => value,
            Bound::Unbounded => panic!("Cannot add unbounded loop"),
        };
        let end = start.max(end);

        match end - start {
            0 => {}
            1 => {
                self.insert(start, start);
            }
            _ => {
                for i in start..end - 1 {
                    self.insert(i, i + 1);
                }
                self.insert(end - 1, start);
            }
        }
    }

    pub fn iter_successors(&self, from: usize) -> impl Iterator<Item = usize> + use<'_> {
        self.0
            .range((from, 0)..=(from, std::usize::MAX))
            .map(|&(_, t)| t)
    }

    pub fn pop_first_outgoing(&mut self, from: usize) -> Option<usize> {
        let to = self.iter_successors(from).next()?;
        self.remove(from, to);
        Some(to)
    }

    pub fn pop_min_outgoing_by<F: Fn(&usize, &usize) -> Ordering>(
        &mut self,
        from: usize,
        cmp: F,
    ) -> Option<usize> {
        let to = self.iter_successors(from).min_by(cmp)?;
        self.remove(from, to);
        Some(to)
    }

    pub fn iter(&self) -> impl Iterator<Item = (usize, usize)> + use<'_> {
        self.0.iter().cloned()
    }

    pub fn len(&self) -> usize {
        self.0.len()
    }
    pub fn contains(&self, from: usize, to: usize) -> bool {
        self.0.contains(&(from, to))
    }
    pub fn retain(&mut self, mut f: impl FnMut(usize, usize) -> bool) {
        self.0.retain(|&(from, to)| f(from, to));
    }
}

impl fmt::Debug for EdgeSet {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "EdgeSet {{")?;
        let mut first_edge = true;
        for &(from, to) in self.0.iter() {
            if !first_edge {
                write!(f, ", ")?;
            } else {
                first_edge = false;
            }
            write!(f, "{:?} -> {:?}", from, to)?;
        }
        write!(f, "}}")?;
        Ok(())
    }
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

fn cmp_avoid(avoid: usize, i: usize, j: usize) -> Ordering {
    (i == avoid).cmp(&(j == avoid))
}

#[derive(Debug, Clone, PartialEq)]
pub enum TriangulationError {
    /// A node with no outgoing edge was detected
    TopologyDeadEnd,

    /// Two coincident points were detected
    CoincidentPoints,

    /// A region with negative area was detected
    /// (e.g. incorrect winding direction)
    NegativeRegion,

    /// A singular point or line segment was detected outside of the polygon
    SingularRegion,

    /// A region or triangle with negative area (or zero area) was detected
    /// (e.g. incorrect winding direction)
    NonPositiveTriangle,

    /// During triangulation, a 1-cycle was detected (1-cycles inside polygons
    /// before they are made monotone are OK)
    Topology1Cycle,

    /// During triangulation, a 2-cycle was detected (2-cycles inside polygons
    /// before they are made monotone are OK)
    Topology2Cycle,

    /// During triangulation, a >3 cycle was detected (graph is not composed
    /// exclusively of 3-cycles),
    TopologyNotTriangle,
}

pub fn partition_into_monotone_components<'a, SPACE: SpaceUv>(
    uv: &[SPACE::Vector],
    mut edges: EdgeSet,
) -> Result<EdgeSet, TriangulationError>
where
    SPACE::Scalar: Debug,
    SPACE::Vector: Debug,
    SPACE::AntiScalar: Debug,
{
    #[derive(Clone)]
    struct MonotoneEdge {
        from_node: usize,
        to_node: usize,
        helper: usize,
        helper_is_merge: bool,
    }
    impl fmt::Debug for MonotoneEdge {
        fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
            write!(f, "{:?} -> {:?}", self.from_node, self.to_node)
        }
    }

    // The state of the sweep line algorithm: a sorted list of all active edges and their helpers.
    // (consider replacing this with a data structure allowing efficient insertion into the middle)
    let mut t = Vec::<MonotoneEdge>::new();

    fn search_t<VECTOR: Copy + Join<VECTOR, Output: Join<VECTOR, Output: Default + PartialOrd>>>(
        t: &Vec<MonotoneEdge>,
        uv: &[VECTOR],
        i: usize,
    ) -> usize {
        let pt = uv[i];

        let comparison = |e: &MonotoneEdge| {
            i != e.from_node
                && i != e.to_node
                && uv[e.from_node].join(uv[e.to_node]).join(pt) > Default::default()
        };

        // Return index of first value that compares >= the given pt
        t.partition_point(comparison)
    }

    fn search_exact_t<
        VECTOR: Copy + Join<VECTOR, Output: Join<VECTOR, Output: Default + PartialOrd>>,
    >(
        t: &Vec<MonotoneEdge>,
        uv: &[VECTOR],
        from: usize,
        to: usize,
    ) -> Option<usize> {
        let mut ix = search_t(t, &uv, to);

        // Within t there may be a range of edges
        // all of which include the given point.
        // ix now points to the first of these,
        // so we can linearly walk forwards until we find one that matches
        // or run out of options.
        while ix < t.len() {
            if t[ix].from_node == from && t[ix].to_node == to {
                return Some(ix);
            }
            ix += 1;
        }
        None
    }

    // Remove the given left edge from the t
    fn pop_t<VECTOR: Copy + Join<VECTOR, Output: Join<VECTOR, Output: Default + PartialOrd>>>(
        t: &mut Vec<MonotoneEdge>,
        uv: &[VECTOR],
        from: usize,
        to: usize,
    ) -> Option<MonotoneEdge> {
        let ix = search_exact_t(t, &uv, from, to)?;
        Some(t.remove(ix))
    }

    // Replace the given left edge from the t
    fn replace_t<
        VECTOR: Copy + Join<VECTOR, Output: Join<VECTOR, Output: Default + PartialOrd>>,
    >(
        t: &mut Vec<MonotoneEdge>,
        uv: &[VECTOR],
        from: usize,
        to: usize,
        new_edge: MonotoneEdge,
    ) -> Option<MonotoneEdge> {
        let ix = search_exact_t(t, &uv, from, to)?;
        Some(replace(&mut t[ix], new_edge))
    }

    fn push_t<VECTOR: Copy + Join<VECTOR, Output: Join<VECTOR, Output: Default + PartialOrd>>>(
        t: &mut Vec<MonotoneEdge>,
        uv: &[VECTOR],
        i: usize,
        edge: MonotoneEdge,
    ) {
        let ix = search_t(t, &uv, i);
        t.insert(ix, edge);
    }

    fn peek_t_mut<
        'a,
        VECTOR: Copy + Join<VECTOR, Output: Join<VECTOR, Output: Default + PartialOrd>>,
    >(
        t: &'a mut Vec<MonotoneEdge>,
        uv: &[VECTOR],
        i: usize,
    ) -> Option<&'a mut MonotoneEdge> {
        let ix = search_t(t, &uv, i).checked_sub(1)?;

        println!(
            "Left of interval is {:?} -> {:?}",
            t[ix].from_node, t[ix].to_node
        );
        Some(&mut t[ix])
    }

    // The order of these is important--
    // it sets the order that the events will be handled
    // if a given vertex has multiple events
    #[derive(Debug, PartialEq, Eq, PartialOrd, Ord)]
    enum EventType {
        End,
        Merge,
        NormalRight,
        NormalLeft,
        Split,
        Start,
        Singular,
    }

    let mut unvisited_edges = edges.clone();

    let event_list = {
        fn categorize_event<
            VECTOR: Sized
                + Copy
                + Dot<VECTOR, Output: Ord>
                + XHat
                + YHat
                + Join<VECTOR, Output: Join<VECTOR, Output: Default + PartialOrd>>,
        >(
            uv_prev: VECTOR,
            uv_cur: VECTOR,
            uv_next: VECTOR,
        ) -> Result<EventType, TriangulationError>
        where
            VECTOR: Debug, // TEMP
        {
            // Figure out the type of vertex i by examining its neighbors
            println!("*** Point {:?} ***", uv_cur);
            println!("Prev is {:?}", uv_prev);
            println!("Next is {:?}", uv_next);

            // Categorize the vertex
            let next_pt_below = match cmp_uv(uv_next, uv_cur) {
                Ordering::Equal => return Err(TriangulationError::CoincidentPoints),
                Ordering::Greater => true,
                Ordering::Less => false,
            };
            let prev_pt_below = match cmp_uv(uv_prev, uv_cur) {
                Ordering::Equal => return Err(TriangulationError::CoincidentPoints),
                Ordering::Greater => true,
                Ordering::Less => false,
            };

            // Calculate the direction we are turning at this vertex
            // to differentiate between start / split vertices and end / merge vertices
            let turn_direction = uv_prev.join(uv_cur).join(uv_next);

            Ok(if next_pt_below && prev_pt_below {
                if turn_direction > Default::default() {
                    EventType::Start
                } else {
                    EventType::Split
                }
            } else if !next_pt_below && !prev_pt_below {
                if turn_direction > Default::default() {
                    EventType::End
                } else {
                    EventType::Merge
                }
            } else {
                if is_on_left_chain(uv_cur, uv_next) {
                    EventType::NormalLeft
                } else {
                    EventType::NormalRight
                }
            })
        }

        let mut event_list = Vec::<(usize, usize, usize, EventType)>::new();
        while let Some((start_node, second_node)) = unvisited_edges.pop() {
            // Walk a single polygon from the graph

            if start_node == second_node {
                println!("Singular point: {:?}", start_node);
                // Special handling for self-loops (singular points)
                event_list.push((start_node, start_node, start_node, EventType::Singular));
                continue;
            }

            println!("START WALK {:?} to {:?}", start_node, second_node);

            let (mut prev_node, mut cur_node) = (start_node, second_node);
            let mut uv_prev = uv[prev_node];
            let mut uv_cur = uv[cur_node];

            loop {
                // Find the next node
                let next_node = unvisited_edges
                    .pop_min_outgoing_by(cur_node, |&i, &j| {
                        cmp_avoid(prev_node, i, j)
                            .then(cmp_turn_angle(uv_prev, uv_cur, uv[i], uv[j]))
                    })
                    .ok_or(TriangulationError::TopologyDeadEnd)?;
                let uv_next = uv[next_node];

                println!("walk from {:?} to {:?}", cur_node, next_node);
                event_list.push((
                    prev_node,
                    cur_node,
                    next_node,
                    categorize_event(uv_prev, uv_cur, uv_next)?,
                ));

                if next_node == start_node {
                    // We didn't actually push the start node because we couldn't determine is
                    // category until we knew its prev
                    event_list.push((
                        cur_node,
                        start_node,
                        second_node,
                        categorize_event(uv_cur, uv_next, uv[second_node])?,
                    ));
                    println!("DONE!");
                    break;
                }

                (prev_node, cur_node) = (cur_node, next_node);
                (uv_prev, uv_cur) = (uv_cur, uv_next);
            }
            // handle loop
        }

        // Sort the event list into sweep-line order
        // (ordering first by -Y coordinate, then by X coordinate, then by event type)
        event_list.sort_by(|(_, i1, _, typ1), (_, i2, _, typ2)| {
            cmp_uv(uv[*i1], uv[*i2]).then(typ1.cmp(typ2))
        });
        event_list
    };

    // Iterate over vertices in sweep-line order
    for (prev, i, next, vertex_type) in event_list.into_iter() {
        println!("* {:?} @ {:?} is a {:?} vertex", i, uv[i], vertex_type);

        // Take action based on the vertex type
        match vertex_type {
            EventType::Start => {
                push_t(
                    &mut t,
                    &uv,
                    i,
                    MonotoneEdge {
                        from_node: i,
                        to_node: next,
                        helper: i,
                        helper_is_merge: false,
                    },
                );
            }
            EventType::Split => {
                {
                    let edge =
                        peek_t_mut(&mut t, &uv, i).ok_or(TriangulationError::NegativeRegion)?;
                    edges.insert(i, edge.helper);
                    edges.insert(edge.helper, i);
                    edge.helper = i;
                    edge.helper_is_merge = false;
                }
                push_t(
                    &mut t,
                    &uv,
                    i,
                    MonotoneEdge {
                        from_node: i,
                        to_node: next,
                        helper: i,
                        helper_is_merge: false,
                    },
                );
            }
            EventType::End => {
                let edge =
                    pop_t(&mut t, &uv, prev, i).expect("End vertex incident edge not found in t");
                println!("Pop {:?}", edge);
                if edge.helper_is_merge {
                    edges.insert(i, edge.helper);
                    edges.insert(edge.helper, i);
                }
            }
            EventType::Merge => {
                let edge =
                    pop_t(&mut t, &uv, prev, i).expect("Merge vertex incident edge not found in t");

                if edge.helper_is_merge {
                    edges.insert(i, edge.helper);
                    edges.insert(edge.helper, i);
                }

                let edge = peek_t_mut(&mut t, &uv, i).ok_or(TriangulationError::NegativeRegion)?;
                if edge.helper_is_merge {
                    edges.insert(i, edge.helper);
                    edges.insert(edge.helper, i);
                }
                edge.helper = i;
                edge.helper_is_merge = true;
            }
            EventType::NormalLeft => {
                // pop exact
                let edge = replace_t(
                    &mut t,
                    &uv,
                    prev,
                    i,
                    MonotoneEdge {
                        from_node: i,
                        to_node: next,
                        helper: i,
                        helper_is_merge: false,
                    },
                )
                .expect("Left vertex incident edge not found in t");
                if edge.helper_is_merge {
                    // Walk the part of the polygon we are cutting off
                    edges.insert(i, edge.helper);
                    edges.insert(edge.helper, i);
                }
            }
            EventType::NormalRight => {
                let edge = peek_t_mut(&mut t, &uv, i).ok_or(TriangulationError::NegativeRegion)?;
                if edge.helper_is_merge {
                    edges.insert(i, edge.helper);
                    edges.insert(edge.helper, i);
                }
                edge.helper = i;
                edge.helper_is_merge = false;
            }
            EventType::Singular => {
                // A singular vertex acts as both a split vertex followed by an immediate merge.
                // Combining those two actions results in this simple result:
                let edge = peek_t_mut(&mut t, &uv, i).ok_or(TriangulationError::SingularRegion)?;
                edges.insert(i, edge.helper);
                edges.insert(edge.helper, i);
                edge.helper = i;
                edge.helper_is_merge = true;
            }
        }
        println!("T is now {:?}", t);
    }

    // Remove any self-loops (they will have been given real new edges by the algorithm)
    edges.retain(|from, to| from != to);

    Ok(edges)
}

pub fn triangulate_monotone_components<SPACE: SpaceUv>(
    uv: &[SPACE::Vector],
    mut edges: EdgeSet,
) -> Result<EdgeSet, TriangulationError>
where
    SPACE::Vector: Debug,     // TEMP
    SPACE::Scalar: Debug,     // TEMP
    SPACE::AntiScalar: Debug, // TEMP
{
    let anti_zero = SPACE::AntiScalar::default();

    // Initialize an "unvisited" set
    let mut unvisited_edges = edges.clone();
    let mut monotone_poly_edges = Vec::<(usize, usize)>::new();
    let mut s = Vec::<(usize, bool)>::new();

    while let Some((start_node, second_node)) = unvisited_edges.pop() {
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

            let next_node = unvisited_edges
                .pop_min_outgoing_by(cur_node, |&i, &j| {
                    cmp_avoid(prev_node, i, j).then(cmp_turn_angle(uv_prev, uv_cur, uv[i], uv[j]))
                })
                .ok_or(TriangulationError::TopologyDeadEnd)?;
            if next_node == cur_node {
                return Err(TriangulationError::SingularRegion);
            }
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
                    s.pop().expect("Stack empty during triangulation");
                println!(
                    "Pop {:?} on the {}",
                    s_top,
                    if is_left { "left" } else { "right" }
                );
                if s_top_is_left != is_left {
                    // Different chains
                    println!("Different chain");
                    println!("Add diagonal {:?} to {:?}", node, s_top);
                    edges.insert(node, s_top);
                    edges.insert(s_top, node);
                    while s.len() > 1 {
                        let (s_top, _) = s.pop().unwrap();

                        // Add an edge to all nodes on the stack,
                        // except the last one
                        println!("Add diagonal {:?} to {:?}", node, s_top);
                        edges.insert(node, s_top);
                        edges.insert(s_top, node);
                    }
                    s.pop().expect("Stack empty during triangulation"); // Discard last stack element
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
                            true => proposed_tri <= anti_zero,
                            false => proposed_tri >= anti_zero,
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
                        edges.insert(node, s_top);
                        edges.insert(s_top, node);
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
                edges.insert(node, s_entry);
                edges.insert(s_entry, node);
            }
        }
    }

    Ok(edges)
}

pub fn edges_to_triangles<SPACE: SpaceUv>(
    uv: &[SPACE::Vector],
    mut edges: EdgeSet,
) -> Result<Vec<[usize; 3]>, TriangulationError> {
    // Initialize an "unvisited" set
    let mut triangles = Vec::<[usize; 3]>::new();

    while let Some((node1, node2)) = edges.pop() {
        // Walk a single monotone polygon from the graph

        println!("*TRI*");
        println!("Node 1: {:?}", node1);
        println!("Node 2: {:?}", node2);
        assert!(node2 != node1, "Graph has a 1-cycle");

        let node3 = {
            let uv_prev = uv[node1];
            let uv_cur = uv[node2];

            let node3 = edges
                .pop_min_outgoing_by(node2, |&i, &j| {
                    let result =
                        cmp_avoid(node1, i, j).then(cmp_turn_angle(uv_prev, uv_cur, uv[i], uv[j]));
                    println!(
                        "CMP: {:?}-{:?}-{:?} {:?} {:?}-{:?}-{:?}",
                        node1, node2, i, result, node1, node2, j
                    );
                    result
                })
                .ok_or(TriangulationError::TopologyDeadEnd)?;

            if node3 == node2 {
                return Err(TriangulationError::Topology1Cycle);
            }
            if node3 == node1 {
                return Err(TriangulationError::Topology2Cycle);
            }
            println!("Node 3: {:?}", node3);
            node3
        };

        {
            let uv_prev = uv[node2];
            let uv_cur = uv[node3];

            let node4 = edges
                .pop_min_outgoing_by(node3, |&i, &j| {
                    let result =
                        cmp_avoid(node2, i, j).then(cmp_turn_angle(uv_prev, uv_cur, uv[i], uv[j]));
                    println!(
                        "CMP: {:?}-{:?}-{:?} {:?} {:?}-{:?}-{:?}",
                        node2, node3, i, result, node2, node3, j
                    );
                    result
                })
                .expect("Bad manifold--no outgoing edge");
            println!("Node 4: {:?}", node4);
            if node4 != node1 {
                return Err(TriangulationError::TopologyNotTriangle);
            }
        }

        if !(uv[node1].join(uv[node2]).join(uv[node3]) > Default::default()) {
            return Err(TriangulationError::NonPositiveTriangle);
        }

        println!("OK");

        triangles.push([node1, node2, node3]);
    }
    Ok(triangles)
}

pub fn triangulate<SPACE: SpaceUv>(
    uv: &[SPACE::Vector],
    edges: EdgeSet,
) -> Result<Vec<[usize; 3]>, TriangulationError>
where
    SPACE::Vector: Debug,     // TEMP
    SPACE::Scalar: Debug,     // TEMP
    SPACE::AntiScalar: Debug, // TEMP
{
    let edges = partition_into_monotone_components::<SPACE>(uv, edges)?;
    let edges = triangulate_monotone_components::<SPACE>(uv, edges)?;
    edges_to_triangles::<SPACE>(uv, edges)
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
    fn test_triangle() {
        let uv = vec![
            Vector::point([R32(0.), R32(0.)]),
            Vector::point([R32(1.), R32(0.)]),
            Vector::point([R32(0.), R32(1.)]),
        ];

        let mut edges = EdgeSet::new();
        edges.insert_loop(0..3);

        let old_edges = edges.clone();
        let edges = partition_into_monotone_components::<Space2D>(&uv, edges).unwrap();

        // Graph should not have changed since a triangle is already monotone
        assert_eq!(old_edges, edges);

        let edges = triangulate_monotone_components::<Space2D>(&uv, edges).unwrap();
        let tri = edges_to_triangles::<Space2D>(&uv, edges).unwrap();
        assert_eq!(tri.len(), 1);
    }

    #[test]
    fn test_square() {
        // Square
        let uv = vec![
            Vector::point([R32(0.), R32(0.)]),
            Vector::point([R32(1.), R32(0.)]),
            Vector::point([R32(1.), R32(1.)]),
            Vector::point([R32(0.), R32(1.)]),
        ];

        let mut edges = EdgeSet::new();
        edges.insert_loop(0..4);

        let old_edges = edges.clone();
        let edges = partition_into_monotone_components::<Space2D>(&uv, edges).unwrap();

        // Graph should not have changed since a triangle is already monotone
        assert_eq!(old_edges, edges);

        let edges = triangulate_monotone_components::<Space2D>(&uv, edges).unwrap();
        let tri = edges_to_triangles::<Space2D>(&uv, edges).unwrap();
        assert_eq!(tri.len(), 2);
    }

    #[test]
    fn test_square_with_redundant_edge_points() {
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

        let mut edges = EdgeSet::new();
        edges.insert_loop(0..8);

        let old_edges = edges.clone();
        let edges = partition_into_monotone_components::<Space2D>(&uv, edges).unwrap();

        // Graph should not have changed since a triangle is already monotone
        assert_eq!(old_edges, edges);

        let edges = triangulate_monotone_components::<Space2D>(&uv, edges).unwrap();
        let tri = edges_to_triangles::<Space2D>(&uv, edges).unwrap();
        assert_eq!(tri.len(), 6);
    }

    #[test]
    fn test_upward_comb() {
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

        let mut edges = EdgeSet::new();
        edges.insert_loop(0..7);

        let edges = partition_into_monotone_components::<Space2D>(&uv, edges).unwrap();

        // Two diagonals should have been added for a total of 4 more graph edges
        assert_eq!(edges.len(), 11);

        let edges = triangulate_monotone_components::<Space2D>(&uv, edges).unwrap();
        let tri = edges_to_triangles::<Space2D>(&uv, edges).unwrap();
        assert_eq!(tri.len(), 5);
    }

    #[test]
    fn test_downward_comb() {
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

        let mut edges = EdgeSet::new();
        edges.insert_loop(0..7);

        let edges = partition_into_monotone_components::<Space2D>(&uv, edges).unwrap();

        // Two diagonals should have been added for a total of 4 more graph edges
        assert_eq!(edges.len(), 11);

        let edges = triangulate_monotone_components::<Space2D>(&uv, edges).unwrap();
        let tri = edges_to_triangles::<Space2D>(&uv, edges).unwrap();
        assert_eq!(tri.len(), 5);
    }

    #[test]
    fn test_square_hole() {
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

        let mut edges = EdgeSet::new();
        edges.insert_loop(0..4);
        edges.insert_loop(4..8);

        let edges = partition_into_monotone_components::<Space2D>(&uv, edges).unwrap();
        // Two diagonals should have been added for a total of 4 more graph edges
        assert_eq!(edges.len(), 12);

        let edges = triangulate_monotone_components::<Space2D>(&uv, edges).unwrap();
        let tri = edges_to_triangles::<Space2D>(&uv, edges).unwrap();
        assert_eq!(tri.len(), 8);
    }

    #[test]
    fn test_multiple_shapes() {
        // Two squares, each with a hole
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

        let mut edges = EdgeSet::new();
        edges.insert_loop(0..4);
        edges.insert_loop(4..8);
        edges.insert_loop(8..12);
        edges.insert_loop(12..16);

        let edges = partition_into_monotone_components::<Space2D>(&uv, edges).unwrap();
        // Four diagonals should have been added for a total of 8 more graph edges
        assert_eq!(edges.len(), 24);

        // Make sure no diagonals go from the first shape (vertices 0-8) to the second shape (vertices 8-16)
        for (from, to) in edges.iter() {
            match from {
                0..8 => assert!((0..8).contains(&to)),
                8..16 => assert!((8..16).contains(&to)),
                _ => panic!(),
            }
        }

        let edges = triangulate_monotone_components::<Space2D>(&uv, edges).unwrap();
        let tri = edges_to_triangles::<Space2D>(&uv, edges).unwrap();
        assert_eq!(tri.len(), 16);
    }

    #[test]
    fn test_line_singularity() {
        // Square with a line inside of it
        let uv = vec![
            // Outside
            Vector::point([R32(0.), R32(0.)]),
            Vector::point([R32(1.), R32(0.)]),
            Vector::point([R32(1.), R32(1.)]),
            Vector::point([R32(0.), R32(1.)]),
            // Two additional points inside it
            Vector::point([R32(0.25), R32(0.5)]),
            Vector::point([R32(0.75), R32(0.5)]),
        ];

        let mut edges = EdgeSet::new();
        edges.insert_loop(0..4);
        // Add a line segment singularity
        edges.insert_loop(4..6);

        let edges = partition_into_monotone_components::<Space2D>(&uv, edges).unwrap();
        // We should have added two diagonals to the singular line
        // to split this into two monotone components
        assert!(edges.contains(0, 5));
        assert!(edges.contains(2, 4));

        let edges = triangulate_monotone_components::<Space2D>(&uv, edges).unwrap();
        let tri = edges_to_triangles::<Space2D>(&uv, edges).unwrap();
        assert_eq!(tri.len(), 6);
    }

    #[test]
    fn test_point_singularity() {
        // Square with a line inside of it
        let uv = vec![
            // Outside
            Vector::point([R32(0.), R32(0.)]),
            Vector::point([R32(1.), R32(0.)]),
            Vector::point([R32(1.), R32(1.)]),
            Vector::point([R32(0.), R32(1.)]),
            // Two additional points inside it
            Vector::point([R32(0.25), R32(0.5)]),
            Vector::point([R32(0.75), R32(0.5)]),
        ];

        let mut edges = EdgeSet::new();
        edges.insert_loop(0..4);
        // Add two point singularities
        edges.insert_loop(4..5);
        edges.insert_loop(5..6);

        let edges = partition_into_monotone_components::<Space2D>(&uv, edges).unwrap();
        // We should have added two diagonals to the singular points
        assert!(edges.contains(0, 5));
        assert!(edges.contains(2, 4));
        // We should have removed the self-loops
        assert!(!edges.contains(4, 4));
        assert!(!edges.contains(5, 5));

        let edges = triangulate_monotone_components::<Space2D>(&uv, edges).unwrap();
        let tri = edges_to_triangles::<Space2D>(&uv, edges).unwrap();
        assert_eq!(tri.len(), 6);
    }

    #[test]
    fn test_starburst() {
        // A square with a 1-dimensional starburst in the middle
        // o--o--o
        // |\ | /|
        // | \|/ |
        // o--o--o
        // | /|\ |
        // |/ | \|
        // o--o--o

        // Square
        let uv = vec![
            // Outside
            Vector::point([R32(0.), R32(0.)]),
            Vector::point([R32(0.5), R32(0.)]),
            Vector::point([R32(1.), R32(0.)]),
            Vector::point([R32(1.), R32(0.5)]),
            Vector::point([R32(1.), R32(1.)]),
            Vector::point([R32(0.5), R32(1.)]),
            Vector::point([R32(0.), R32(1.)]),
            Vector::point([R32(0.), R32(0.5)]),
            // Center
            Vector::point([R32(0.5), R32(0.5)]),
        ];

        let mut edges = EdgeSet::new();
        edges.insert_loop(0..8);
        edges.insert(0, 8);
        edges.insert(1, 8);
        edges.insert(2, 8);
        edges.insert(3, 8);
        edges.insert(4, 8);
        edges.insert(5, 8);
        edges.insert(6, 8);
        edges.insert(7, 8);
        edges.insert(8, 0);
        edges.insert(8, 1);
        edges.insert(8, 2);
        edges.insert(8, 3);
        edges.insert(8, 4);
        edges.insert(8, 5);
        edges.insert(8, 6);
        edges.insert(8, 7);

        let edges = partition_into_monotone_components::<Space2D>(&uv, edges).unwrap();
        // Two diagonals should have been added for a total of 4 more graph edges
        assert_eq!(edges.len(), 24);

        let edges = triangulate_monotone_components::<Space2D>(&uv, edges).unwrap();
        let tri = edges_to_triangles::<Space2D>(&uv, edges).unwrap();
        assert_eq!(tri.len(), 8);
    }

    #[test]
    fn test_make_monotone_errors() {
        // Square
        let uv = vec![
            Vector::point([R32(0.), R32(0.)]),
            Vector::point([R32(1.), R32(0.)]),
            Vector::point([R32(1.), R32(1.)]),
            Vector::point([R32(0.), R32(1.)]),
        ];

        {
            // Branching
            let mut edges = EdgeSet::new();
            edges.insert_loop(0..4);
            edges.insert(0, 2);

            let err = partition_into_monotone_components::<Space2D>(&uv, edges).unwrap_err();
            assert_eq!(err, TriangulationError::TopologyDeadEnd);
        }
        {
            // Dead end
            let mut edges = EdgeSet::new();
            edges.insert(0, 1);
            edges.insert(1, 2);
            edges.insert(2, 3);

            let err = partition_into_monotone_components::<Space2D>(&uv, edges).unwrap_err();
            assert_eq!(err, TriangulationError::TopologyDeadEnd);
        }
        {
            // Wound backwards
            let mut edges = EdgeSet::new();
            edges.insert(3, 2);
            edges.insert(2, 1);
            edges.insert(1, 0);
            edges.insert(0, 3);

            let err = partition_into_monotone_components::<Space2D>(&uv, edges).unwrap_err();
            assert_eq!(err, TriangulationError::NegativeRegion);
        }
        {
            // Self-intersecting
            let mut edges = EdgeSet::new();
            edges.insert(1, 0);
            edges.insert(0, 2);
            edges.insert(2, 3);
            edges.insert(3, 1);

            let err = partition_into_monotone_components::<Space2D>(&uv, edges).unwrap_err();
            assert_eq!(err, TriangulationError::NegativeRegion);
        }
        {
            // Triangle with redundant point
            let uv = vec![
                Vector::point([R32(0.), R32(0.)]),
                Vector::point([R32(1.), R32(0.)]),
                Vector::point([R32(1.), R32(1.)]),
                Vector::point([R32(0.), R32(0.)]),
            ];

            let mut edges = EdgeSet::new();
            edges.insert_loop(0..4);

            let err = partition_into_monotone_components::<Space2D>(&uv, edges).unwrap_err();
            assert_eq!(err, TriangulationError::CoincidentPoints);
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
