use std::{fmt::Display, ops::Range};

use crate::data::ranges::InsertionIdx;

/// A helper struct for incrementally combining ranges/indices into a
/// [`String`].
#[derive(Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug)]
pub struct Coords(String);

impl Coords {
    /// Creates a new [`Coords`] such that the underlying [`String`] has the
    /// specified `capacity`.
    pub(crate) fn with_capacity(capacity: usize) -> Self {
        Self(String::with_capacity(capacity))
    }

    /// PUshes a range to the coordinates.
    pub(crate) fn push_range(&mut self, range: &Range<usize>) {
        if !self.0.is_empty() {
            self.0.push(';');
        }

        let mut buff = core::fmt::NumBuffer::new();

        // 0-based half-open to 1-based inclusive
        self.0.push_str((range.start + 1).format_into(&mut buff));
        self.0.push_str("..");
        self.0.push_str(range.end.format_into(&mut buff));
    }

    // TODO: Rename?

    /// Pushes an insertion index to the coordinates.
    pub(crate) fn push_upstream(&mut self, index: InsertionIdx) {
        if !self.0.is_empty() {
            self.0.push(';');
        }

        let mut buff = core::fmt::NumBuffer::new();

        self.0.push_str(index.left_pos().format_into(&mut buff));
    }
}

impl Display for Coords {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.0.fmt(f)
    }
}

impl AsRef<str> for Coords {
    fn as_ref(&self) -> &str {
        &self.0
    }
}
