use ngeom::algebraic_ops::*;
use ngeom::ops::*;
use ngeom::scalar::*;
use std::cmp::Ordering;
use std::fmt;
use std::fmt::Debug;
use std::mem::replace;
use std::ops::Mul;

use crate::graph::{EdgeSet, NodeIndex};

fn u<VECTOR: Sized + Dot<VECTOR> + XHat>(pt: VECTOR) -> <VECTOR as Dot<VECTOR>>::Output {
    pt.dot(VECTOR::x_hat())
}
fn v<VECTOR: Sized + Dot<VECTOR> + YHat>(pt: VECTOR) -> <VECTOR as Dot<VECTOR>>::Output {
    pt.dot(VECTOR::y_hat())
}
fn partial_cmp_uv<VECTOR: Sized + Copy + Dot<VECTOR, Output: PartialOrd> + XHat + YHat>(
    uv1: VECTOR,
    uv2: VECTOR,
) -> Option<Ordering> {
    Some(
        v(uv1)
            .partial_cmp(&v(uv2))?
            .reverse()
            .then(u(uv1).partial_cmp(&u(uv2))?),
    )
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

fn partial_cmp_turn_angle<
    SCALAR: Copy + Default + Mul<SCALAR, Output = SCALAR>,
    VECTOR: Copy + Join<VECTOR, Output = BIVECTOR>,
    BIVECTOR: Copy
        + Join<VECTOR, Output = ANTISCALAR>
        + SubsetOrthogonalTo<BIVECTOR, Output = SCALAR>
        + WeightNorm<Output = ANTISCALAR>,
    ANTISCALAR: Copy + Default + Mul<SCALAR, Output = ANTISCALAR> + PartialOrd,
>(
    pt1: VECTOR,
    pt2: VECTOR,
    pt3a: VECTOR,
    pt3b: VECTOR,
) -> Option<Ordering> {
    let l1 = pt1.join(pt2);
    let l2a = pt2.join(pt3a);
    let l2b = pt2.join(pt3b);

    let a_w = l2a.weight_norm();
    let b_w = l2b.weight_norm();
    let a_sin = l1.join(pt3a);
    let b_sin = l1.join(pt3b);
    let a_cos = l1.subset_orthogonal_to(l2a);
    let b_cos = l1.subset_orthogonal_to(l2b);

    // For a valid comparison, we must normalize by the outgoing segment weight
    // which may be different between options A and B.
    // We can do this without division by moving the weight to the opposite side of the equation,
    // since both weights are positive.

    if a_sin >= Default::default() && b_sin >= Default::default() {
        // Both angles are not reflex,
        // so compare their cosines
        (b_w * a_cos).partial_cmp(&(a_w * b_cos))
    } else if a_sin < Default::default() && b_sin < Default::default() {
        // Both angles are reflex
        // so compare their cosines
        (a_w * b_cos).partial_cmp(&(b_w * a_cos))
    } else {
        // One angle is reflex and one is normal.
        // That means one sine is positive and one is non-negative.
        // This trivial comparison will return the appropriate result
        b_sin.partial_cmp(&a_sin)
    }
}

fn cmp_avoid(avoid: NodeIndex, i: NodeIndex, j: NodeIndex) -> Ordering {
    (i == avoid).cmp(&(j == avoid))
}

fn walk_sharpest_angle<
    SCALAR: Copy + Default + Mul<SCALAR, Output = SCALAR>,
    VECTOR: Copy + Join<VECTOR, Output = BIVECTOR>,
    BIVECTOR: Copy
        + Join<VECTOR, Output = ANTISCALAR>
        + SubsetOrthogonalTo<BIVECTOR, Output = SCALAR>
        + WeightNorm<Output = ANTISCALAR>,
    ANTISCALAR: Copy + Default + Mul<SCALAR, Output = ANTISCALAR> + PartialOrd,
>(
    uv: &[VECTOR],
    unvisited_edges: &mut EdgeSet,
    prev_node: NodeIndex,
    cur_node: NodeIndex,
) -> Result<NodeIndex, TriangulationError> {
    let uv_prev = uv[prev_node];
    let uv_cur = uv[cur_node];

    let mut bad_cmp = false;
    let next = unvisited_edges
        .pop_min_outgoing_by(cur_node, |&i, &j| {
            cmp_avoid(prev_node, i, j).then_with(|| {
                partial_cmp_turn_angle(uv_prev, uv_cur, uv[i], uv[j]).unwrap_or_else(|| {
                    bad_cmp = true;
                    Ordering::Equal
                })
            })
        })
        .ok_or(TriangulationError::TopologyDeadEnd);

    if bad_cmp {
        return Err(TriangulationError::IncomparablePoints);
    }

    next
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

    /// During triangulation, two points were found whose coordinates could not be compared
    /// (e.g. were NAN)
    IncomparablePoints,
}

pub fn partition_into_monotone_components<
    'a,
    SCALAR: Copy + Ring + PartialOrd,
    VECTOR: Copy + Join<VECTOR, Output = BIVECTOR> + Dot<VECTOR, Output = SCALAR> + XHat + YHat,
    BIVECTOR: Copy
        + Join<VECTOR, Output = ANTISCALAR>
        + SubsetOrthogonalTo<BIVECTOR, Output = SCALAR>
        + WeightNorm<Output = ANTISCALAR>,
    ANTISCALAR: Copy + Default + PartialOrd + Mul<SCALAR, Output = ANTISCALAR>,
>(
    uv: &[VECTOR],
    mut edges: EdgeSet,
) -> Result<EdgeSet, TriangulationError>
where
    //Scalar: Debug,
    VECTOR: Debug,
    //AntiScalar: Debug,
{
    #[derive(Clone)]
    struct MonotoneEdge {
        from_node: NodeIndex,
        to_node: NodeIndex,
        helper: NodeIndex,
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
        i: NodeIndex,
    ) -> NodeIndex {
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
        from: NodeIndex,
        to: NodeIndex,
    ) -> Option<NodeIndex> {
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
        from: NodeIndex,
        to: NodeIndex,
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
        from: NodeIndex,
        to: NodeIndex,
        new_edge: MonotoneEdge,
    ) -> Option<MonotoneEdge> {
        let ix = search_exact_t(t, &uv, from, to)?;
        Some(replace(&mut t[ix], new_edge))
    }

    fn push_t<VECTOR: Copy + Join<VECTOR, Output: Join<VECTOR, Output: Default + PartialOrd>>>(
        t: &mut Vec<MonotoneEdge>,
        uv: &[VECTOR],
        i: NodeIndex,
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
        i: NodeIndex,
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
                + Dot<VECTOR, Output: PartialOrd>
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
            let next_pt_below = match partial_cmp_uv(uv_next, uv_cur)
                .ok_or(TriangulationError::IncomparablePoints)?
            {
                Ordering::Equal => return Err(TriangulationError::CoincidentPoints),
                Ordering::Greater => true,
                Ordering::Less => false,
            };
            let prev_pt_below = match partial_cmp_uv(uv_prev, uv_cur)
                .ok_or(TriangulationError::IncomparablePoints)?
            {
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

        let mut event_list = Vec::<(NodeIndex, NodeIndex, NodeIndex, EventType)>::new();
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
                let next_node =
                    walk_sharpest_angle(&uv, &mut unvisited_edges, prev_node, cur_node)?;
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
        let mut bad_cmp = false;
        event_list.sort_by(|(_, i1, _, typ1), (_, i2, _, typ2)| {
            partial_cmp_uv(uv[*i1], uv[*i2])
                .unwrap_or_else(|| {
                    bad_cmp = true;
                    Ordering::Equal
                })
                .then(typ1.cmp(typ2))
        });
        if bad_cmp {
            return Err(TriangulationError::IncomparablePoints);
        }
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

pub fn triangulate_monotone_components<
    SCALAR: Copy + Ring + PartialOrd,
    VECTOR: Copy + Join<VECTOR, Output = BIVECTOR> + Dot<VECTOR, Output = SCALAR> + XHat + YHat,
    BIVECTOR: Copy
        + Join<VECTOR, Output = ANTISCALAR>
        + SubsetOrthogonalTo<BIVECTOR, Output = SCALAR>
        + WeightNorm<Output = ANTISCALAR>,
    ANTISCALAR: Copy + Default + PartialOrd + Mul<SCALAR, Output = ANTISCALAR>,
>(
    uv: &[VECTOR],
    mut edges: EdgeSet,
) -> Result<EdgeSet, TriangulationError>
where
    SCALAR: Debug,     // TEMP
    VECTOR: Debug,     // TEMP
    BIVECTOR: Debug,   // TEMP
    ANTISCALAR: Debug, // TEMP
{
    // Initialize an "unvisited" set
    let mut unvisited_edges = edges.clone();
    let mut monotone_poly_edges = Vec::<(NodeIndex, NodeIndex)>::new();
    let mut s = Vec::<(NodeIndex, bool)>::new();

    while let Some((start_node, second_node)) = unvisited_edges.pop() {
        // Walk a single monotone polygon from the graph
        monotone_poly_edges.clear(); // Prepare for a fresh polygon
        monotone_poly_edges.push((start_node, second_node));

        println!("START WALK {:?} to {:?}", start_node, second_node);

        let (mut prev_node, mut cur_node) = (start_node, second_node);

        loop {
            println!("walk from {:?} to {:?}", prev_node, cur_node);

            // Find the next node
            let next_node = walk_sharpest_angle(&uv, &mut unvisited_edges, prev_node, cur_node)?;
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
            let mut bad_cmp = false;
            monotone_poly_edges.sort_by(|&(i, _), &(j, _)| {
                partial_cmp_uv(uv[i], uv[j]).unwrap_or_else(|| {
                    bad_cmp = true;
                    Ordering::Equal
                })
            });
            if bad_cmp {
                return Err(TriangulationError::IncomparablePoints);
            }

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
                            true => proposed_tri <= Default::default(),
                            false => proposed_tri >= Default::default(),
                        } {
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

pub fn edges_to_triangles<
    SCALAR: Copy + Ring + PartialOrd,
    VECTOR: Copy + Join<VECTOR, Output = BIVECTOR> + Dot<VECTOR, Output = SCALAR> + XHat + YHat,
    BIVECTOR: Copy
        + Join<VECTOR, Output = ANTISCALAR>
        + SubsetOrthogonalTo<BIVECTOR, Output = SCALAR>
        + WeightNorm<Output = ANTISCALAR>,
    ANTISCALAR: Copy + Default + PartialOrd + Mul<SCALAR, Output = ANTISCALAR>,
>(
    uv: &[VECTOR],
    mut edges: EdgeSet,
) -> Result<Vec<[NodeIndex; 3]>, TriangulationError> {
    // Initialize an "unvisited" set
    let mut triangles = Vec::<[NodeIndex; 3]>::new();

    while let Some((node1, node2)) = edges.pop() {
        // Walk a single monotone polygon from the graph

        println!("*TRI*");
        println!("Node 1: {:?}", node1);
        println!("Node 2: {:?}", node2);
        assert!(node2 != node1, "Graph has a 1-cycle");

        let node3 = {
            let node3 = walk_sharpest_angle(&uv, &mut edges, node1, node2)?;

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
            let node4 = walk_sharpest_angle(&uv, &mut edges, node2, node3)?;
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

pub fn triangulate<
    SCALAR: Copy + Ring + PartialOrd,
    VECTOR: Copy + Join<VECTOR, Output = BIVECTOR> + Dot<VECTOR, Output = SCALAR> + XHat + YHat,
    BIVECTOR: Copy
        + Join<VECTOR, Output = ANTISCALAR>
        + SubsetOrthogonalTo<BIVECTOR, Output = SCALAR>
        + WeightNorm<Output = ANTISCALAR>,
    ANTISCALAR: Copy + Default + PartialOrd + Mul<SCALAR, Output = ANTISCALAR>,
>(
    uv: &[VECTOR],
    edges: EdgeSet,
) -> Result<Vec<[NodeIndex; 3]>, TriangulationError>
where
    SCALAR: Debug,     // TEMP
    VECTOR: Debug,     // TEMP
    BIVECTOR: Debug,   // TEMP
    ANTISCALAR: Debug, // TEMP
{
    let edges = partition_into_monotone_components(uv, edges)?;
    let edges = triangulate_monotone_components(uv, edges)?;
    edges_to_triangles(uv, edges)
}

#[cfg(test)]
mod tests {
    use ngeom::re2::Vector;

    use super::*;

    #[test]
    fn test_triangle() {
        let uv = vec![
            Vector::point([0., 0.]),
            Vector::point([1., 0.]),
            Vector::point([0., 1.]),
        ];

        let mut edges = EdgeSet::new();
        edges.insert_loop(0..3);

        let old_edges = edges.clone();
        let edges = partition_into_monotone_components(&uv, edges).unwrap();

        // Graph should not have changed since a triangle is already monotone
        assert_eq!(old_edges, edges);

        let edges = triangulate_monotone_components(&uv, edges).unwrap();
        let tri = edges_to_triangles(&uv, edges).unwrap();
        assert_eq!(tri.len(), 1);
    }

    #[test]
    fn test_square() {
        // Square
        let uv = vec![
            Vector::point([0., 0.]),
            Vector::point([1., 0.]),
            Vector::point([1., 1.]),
            Vector::point([0., 1.]),
        ];

        let mut edges = EdgeSet::new();
        edges.insert_loop(0..4);

        let old_edges = edges.clone();
        let edges = partition_into_monotone_components(&uv, edges).unwrap();

        // Graph should not have changed since a triangle is already monotone
        assert_eq!(old_edges, edges);

        let edges = triangulate_monotone_components(&uv, edges).unwrap();
        let tri = edges_to_triangles(&uv, edges).unwrap();
        assert_eq!(tri.len(), 2);
    }

    #[test]
    fn test_square_with_redundant_edge_points() {
        // Square with redundant edge points
        let uv = vec![
            Vector::point([0., 0.]),
            Vector::point([0.5, 0.]),
            Vector::point([1., 0.]),
            Vector::point([1., 0.5]),
            Vector::point([1., 1.]),
            Vector::point([0.5, 1.]),
            Vector::point([0., 1.]),
            Vector::point([0., 0.5]),
        ];

        let mut edges = EdgeSet::new();
        edges.insert_loop(0..8);

        let old_edges = edges.clone();
        let edges = partition_into_monotone_components(&uv, edges).unwrap();

        // Graph should not have changed since a triangle is already monotone
        assert_eq!(old_edges, edges);

        let edges = triangulate_monotone_components(&uv, edges).unwrap();
        let tri = edges_to_triangles(&uv, edges).unwrap();
        assert_eq!(tri.len(), 6);
    }

    #[test]
    fn test_upward_comb() {
        // Three-pointed comb like this: /\/\/\
        let uv = vec![
            Vector::point([6., 0.]),
            Vector::point([5., 1.]),
            Vector::point([4., 0.5]),
            Vector::point([3., 1.]),
            Vector::point([2., 0.5]),
            Vector::point([1., 1.]),
            Vector::point([0., 0.]),
        ];

        let mut edges = EdgeSet::new();
        edges.insert_loop(0..7);

        let edges = partition_into_monotone_components(&uv, edges).unwrap();

        // Two diagonals should have been added for a total of 4 more graph edges
        assert_eq!(edges.len(), 11);

        let edges = triangulate_monotone_components(&uv, edges).unwrap();
        let tri = edges_to_triangles(&uv, edges).unwrap();
        assert_eq!(tri.len(), 5);
    }

    #[test]
    fn test_downward_comb() {
        // Upside-down three-pointed comb like this: \/\/\/
        let uv = vec![
            Vector::point([0., 0.]),
            Vector::point([1., -1.]),
            Vector::point([2., -0.5]),
            Vector::point([3., -1.]),
            Vector::point([4., -0.5]),
            Vector::point([5., -1.]),
            Vector::point([6., 0.]),
        ];

        let mut edges = EdgeSet::new();
        edges.insert_loop(0..7);

        let edges = partition_into_monotone_components(&uv, edges).unwrap();

        // Two diagonals should have been added for a total of 4 more graph edges
        assert_eq!(edges.len(), 11);

        let edges = triangulate_monotone_components(&uv, edges).unwrap();
        let tri = edges_to_triangles(&uv, edges).unwrap();
        assert_eq!(tri.len(), 5);
    }

    #[test]
    fn test_square_hole() {
        // Square
        let uv = vec![
            // Outside
            Vector::point([0., 0.]),
            Vector::point([1., 0.]),
            Vector::point([1., 1.]),
            Vector::point([0., 1.]),
            // Hole
            Vector::point([0.25, 0.75]),
            Vector::point([0.75, 0.75]),
            Vector::point([0.75, 0.25]),
            Vector::point([0.25, 0.25]),
        ];

        let mut edges = EdgeSet::new();
        edges.insert_loop(0..4);
        edges.insert_loop(4..8);

        let edges = partition_into_monotone_components(&uv, edges).unwrap();
        // Two diagonals should have been added for a total of 4 more graph edges
        assert_eq!(edges.len(), 12);

        let edges = triangulate_monotone_components(&uv, edges).unwrap();
        let tri = edges_to_triangles(&uv, edges).unwrap();
        assert_eq!(tri.len(), 8);
    }

    #[test]
    fn test_multiple_shapes() {
        // Two squares, each with a hole
        let uv = vec![
            // Outside
            Vector::point([0., 0.]),
            Vector::point([1., 0.]),
            Vector::point([1., 1.]),
            Vector::point([0., 1.]),
            // Hole
            Vector::point([0.25, 0.75]),
            Vector::point([0.75, 0.75]),
            Vector::point([0.75, 0.25]),
            Vector::point([0.25, 0.25]),
            // Second square outside
            Vector::point([2., 0.]),
            Vector::point([3., 0.]),
            Vector::point([3., 1.]),
            Vector::point([2., 1.]),
            // Second square hole
            Vector::point([2.25, 0.75]),
            Vector::point([2.75, 0.75]),
            Vector::point([2.75, 0.25]),
            Vector::point([2.25, 0.25]),
        ];

        let mut edges = EdgeSet::new();
        edges.insert_loop(0..4);
        edges.insert_loop(4..8);
        edges.insert_loop(8..12);
        edges.insert_loop(12..16);

        let edges = partition_into_monotone_components(&uv, edges).unwrap();
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

        let edges = triangulate_monotone_components(&uv, edges).unwrap();
        let tri = edges_to_triangles(&uv, edges).unwrap();
        assert_eq!(tri.len(), 16);
    }

    #[test]
    fn test_line_singularity() {
        // Square with a line inside of it
        let uv = vec![
            // Outside
            Vector::point([0., 0.]),
            Vector::point([1., 0.]),
            Vector::point([1., 1.]),
            Vector::point([0., 1.]),
            // Two additional points inside it
            Vector::point([0.25, 0.5]),
            Vector::point([0.75, 0.5]),
        ];

        let mut edges = EdgeSet::new();
        edges.insert_loop(0..4);
        // Add a line segment singularity
        edges.insert_loop(4..6);

        let edges = partition_into_monotone_components(&uv, edges).unwrap();
        // We should have added two diagonals to the singular line
        // to split this into two monotone components
        assert!(edges.contains(0, 5));
        assert!(edges.contains(2, 4));

        let edges = triangulate_monotone_components(&uv, edges).unwrap();
        let tri = edges_to_triangles(&uv, edges).unwrap();
        assert_eq!(tri.len(), 6);
    }

    #[test]
    fn test_point_singularity() {
        // Square with a line inside of it
        let uv = vec![
            // Outside
            Vector::point([0., 0.]),
            Vector::point([1., 0.]),
            Vector::point([1., 1.]),
            Vector::point([0., 1.]),
            // Two additional points inside it
            Vector::point([0.25, 0.5]),
            Vector::point([0.75, 0.5]),
        ];

        let mut edges = EdgeSet::new();
        edges.insert_loop(0..4);
        // Add two point singularities
        edges.insert_loop(4..5);
        edges.insert_loop(5..6);

        let edges = partition_into_monotone_components(&uv, edges).unwrap();
        // We should have added two diagonals to the singular points
        assert!(edges.contains(0, 5));
        assert!(edges.contains(2, 4));
        // We should have removed the self-loops
        assert!(!edges.contains(4, 4));
        assert!(!edges.contains(5, 5));

        let edges = triangulate_monotone_components(&uv, edges).unwrap();
        let tri = edges_to_triangles(&uv, edges).unwrap();
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
            Vector::point([0., 0.]),
            Vector::point([0.5, 0.]),
            Vector::point([1., 0.]),
            Vector::point([1., 0.5]),
            Vector::point([1., 1.]),
            Vector::point([0.5, 1.]),
            Vector::point([0., 1.]),
            Vector::point([0., 0.5]),
            // Center
            Vector::point([0.5, 0.5]),
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

        let edges = partition_into_monotone_components(&uv, edges).unwrap();
        // Two diagonals should have been added for a total of 4 more graph edges
        assert_eq!(edges.len(), 24);

        let edges = triangulate_monotone_components(&uv, edges).unwrap();
        let tri = edges_to_triangles(&uv, edges).unwrap();
        assert_eq!(tri.len(), 8);
    }

    #[test]
    fn test_make_monotone_errors() {
        // Square
        let uv = vec![
            Vector::point([0., 0.]),
            Vector::point([1., 0.]),
            Vector::point([1., 1.]),
            Vector::point([0., 1.]),
        ];

        {
            // Branching
            let mut edges = EdgeSet::new();
            edges.insert_loop(0..4);
            edges.insert(0, 2);

            let err = partition_into_monotone_components(&uv, edges).unwrap_err();
            assert_eq!(err, TriangulationError::TopologyDeadEnd);
        }
        {
            // Dead end
            let mut edges = EdgeSet::new();
            edges.insert(0, 1);
            edges.insert(1, 2);
            edges.insert(2, 3);

            let err = partition_into_monotone_components(&uv, edges).unwrap_err();
            assert_eq!(err, TriangulationError::TopologyDeadEnd);
        }
        {
            // Wound backwards
            let mut edges = EdgeSet::new();
            edges.insert(3, 2);
            edges.insert(2, 1);
            edges.insert(1, 0);
            edges.insert(0, 3);

            let err = partition_into_monotone_components(&uv, edges).unwrap_err();
            assert_eq!(err, TriangulationError::NegativeRegion);
        }
        {
            // Self-intersecting
            let mut edges = EdgeSet::new();
            edges.insert(1, 0);
            edges.insert(0, 2);
            edges.insert(2, 3);
            edges.insert(3, 1);

            let err = partition_into_monotone_components(&uv, edges).unwrap_err();
            assert_eq!(err, TriangulationError::NegativeRegion);
        }
        {
            // Triangle with redundant point
            let uv = vec![
                Vector::point([0., 0.]),
                Vector::point([1., 0.]),
                Vector::point([1., 1.]),
                Vector::point([0., 0.]),
            ];

            let mut edges = EdgeSet::new();
            edges.insert_loop(0..4);

            let err = partition_into_monotone_components(&uv, edges).unwrap_err();
            assert_eq!(err, TriangulationError::CoincidentPoints);
        }
    }
}
