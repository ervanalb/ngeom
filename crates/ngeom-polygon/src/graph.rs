use std::cmp::Ordering;
use std::collections::BTreeSet;
use std::fmt;
use std::iter::Iterator;
use std::ops::{Bound, RangeBounds};

pub type NodeIndex = usize;

#[derive(Clone, Default, PartialEq, Eq)]
pub struct EdgeSet(BTreeSet<(NodeIndex, NodeIndex)>);

impl EdgeSet {
    pub fn new() -> Self {
        EdgeSet(BTreeSet::new())
    }
    pub fn insert(&mut self, from: NodeIndex, to: NodeIndex) -> bool {
        self.0.insert((from, to))
    }
    pub fn remove(&mut self, from: NodeIndex, to: NodeIndex) -> bool {
        self.0.remove(&(from, to))
    }
    pub fn pop(&mut self) -> Option<(NodeIndex, NodeIndex)> {
        self.0.pop_last()
    }
    pub fn insert_loop(&mut self, range: impl RangeBounds<NodeIndex>) {
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

    pub fn iter_successors(&self, from: NodeIndex) -> impl Iterator<Item = NodeIndex> + use<'_> {
        self.0.range((from, 0)..(from + 1, 0)).map(|&(_, t)| t)
    }

    pub fn pop_first_outgoing(&mut self, from: NodeIndex) -> Option<NodeIndex> {
        let to = self.iter_successors(from).next()?;
        self.remove(from, to);
        Some(to)
    }

    pub fn pop_min_outgoing_by(
        &mut self,
        from: NodeIndex,
        cmp: impl FnMut(&NodeIndex, &NodeIndex) -> Ordering,
    ) -> Option<NodeIndex> {
        let to = self.iter_successors(from).min_by(cmp)?;

        self.remove(from, to);
        Some(to)
    }

    pub fn iter(&self) -> impl Iterator<Item = (NodeIndex, NodeIndex)> + use<'_> {
        self.0.iter().cloned()
    }

    pub fn len(&self) -> NodeIndex {
        self.0.len()
    }
    pub fn contains(&self, from: NodeIndex, to: NodeIndex) -> bool {
        self.0.contains(&(from, to))
    }
    pub fn retain(&mut self, mut f: impl FnMut(NodeIndex, NodeIndex) -> bool) {
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
