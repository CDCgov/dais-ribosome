//! Structs and enums for representing product and genome alignments,
//! specifically the ranges where matches, insertions, and deletions apply.

use crate::tsv::{Nullable, NullableValue};
use std::{
    cmp::Ordering,
    fmt::{self, Display, Formatter},
    marker::PhantomData,
    ops::Range,
};
use zoe::{
    alignment::Alignment,
    data::{cigar::Ciglet, err::ResultWithErrorContext},
};

/// A helper struct to avoid confusion when storing the 0-based index of an
/// insertion within a reference or coding sequence.
///
/// Confusion can occur because we might be using the index before the insertion
/// or the index after the insertion (since insertions happen _between_
/// indices). This struct abstracts this logic, offering better correctness.
#[repr(transparent)]
#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug)]
pub struct InsertionIdx {
    /// The index that occurs after the insertion. Equivalently, the insertion
    /// happens before this index.
    right: usize,
}

impl InsertionIdx {
    /// Creates a new [`InsertionIdx`] representing an insertion between index
    /// `left` and `right=left+1`.
    #[allow(dead_code)]
    pub fn from_left_idx(left: usize) -> Self {
        Self { right: left + 1 }
    }

    /// Creates a new [`InsertionIdx`] representing an insertion between index
    /// `left=right-1` and `right`.
    pub fn from_right_idx(right: usize) -> Self {
        Self { right }
    }

    /// Creates a new [`InsertionIdx`] representing an insertion between the
    /// 1-based positions `left` and `right=left+1`.
    pub fn from_left_pos(left: usize) -> Self {
        // 1-based left is equivalent to 0-based right
        Self { right: left }
    }

    /// Creates a new [`InsertionIdx`] representing an insertion between the
    /// 1-based positions `left=right-1` and `right`.
    ///
    /// ## Panics
    ///
    /// This will panic if `right` is 0.
    pub fn from_right_pos(right: usize) -> Self {
        // Convert 1-based to 0-based
        Self { right: right - 1 }
    }

    /// The 0-based index in the reference sequence before the insertion.
    ///
    /// ## Panics
    ///
    /// This will panic if the insertion occurs at the beginning of the sequence
    /// (right index 0).
    #[allow(dead_code)]
    pub(crate) fn left(self) -> usize {
        self.right - 1
    }

    /// The 0-based index in the reference sequence after the insertion.
    pub(crate) fn right(self) -> usize {
        self.right
    }

    /// A mutable reference to the 0-based index in the reference sequence after
    /// the insertion.
    pub(crate) fn right_mut(&mut self) -> &mut usize {
        &mut self.right
    }

    /// Converts a [`InsertionIdx`] in a nucleotides sequences to the
    /// corresponding [`InsertionIdx`] in the amino acid sequence.
    ///
    /// This should be called only on [`InsertionIdx`] within a nucleotides
    /// sequence (not an amino acids sequence).
    ///
    /// If the insertion appears within a codon rather than between codons, then
    /// it is said that the insertion occurs before that codon's index. In other
    /// words, the right amino acid index of the insertion is the codon's index.
    pub(crate) fn to_aa_idx(self) -> Self {
        // nt_right = 0  ->  insertion before first codon  ->  aa_right = 0
        // nt_right = 1  ->  insertion in first codon      ->  aa_right = 0
        // nt_right = 2  ->  insertion in first codon      ->  aa_right = 0
        // nt_right = 3  ->  insertion before second codon ->  aa_right = 1
        Self { right: self.right / 3 }
    }

    /// Computes the "codon shift" of the insertion.
    ///
    /// See [`ComputedInsertion::codon_shift`].
    ///
    /// ## Validity
    ///
    /// This should be called only on [`InsertionIdx`] within a nucleotides
    /// sequence (not an amino acids sequence).
    ///
    /// [`ComputedInsertion::codon_shift`]:
    ///     crate::outputs::ComputedInsertion::codon_shift
    pub(crate) fn codon_shift(self) -> usize {
        // right = 0  ->  in frame
        // right = 1  ->  after first base
        // right = 2  ->  after second base
        // right = 3  ->  in frame
        self.right % 3
    }

    /// The 1-based position in the reference sequence before the insertion.
    pub(crate) fn left_pos(self) -> usize {
        // left_pos = left_idx + 1 = right_idx
        self.right
    }

    /// The 1-based position in the reference sequence after the insertion.
    pub(crate) fn right_pos(self) -> usize {
        // right_pos = right_idx + 1
        self.right + 1
    }

    /// Returns whether the insertion occurs at the very start of the sequence.
    pub(crate) fn at_start(self) -> bool {
        // right = 0  ->  insertion before first base
        // right = 1  ->  insertion after first base
        self.right == 0
    }
}

/// The range/index within the query and reference where a contiguous block of
/// matches, deletions, or insertions occur.
///
/// This is a helper type for protein annotation and coordinate manipulation.
#[derive(Clone, Debug)]
pub enum StateRange {
    M(MatchRange),
    D(DeletionRange),
    I(InsertionRange),
}

impl StateRange {
    /// Returns the length of the state in either the query or reference
    /// coordinates, whichever is applicable.
    pub fn len(&self) -> usize {
        match self {
            StateRange::M(range) => range.len(),
            StateRange::D(range) => range.len(),
            StateRange::I(range) => range.len(),
        }
    }

    /// Returns whether the state is empty (zero length) in either the query or
    /// reference coordinates, whichever is applicable.
    pub fn is_empty(&self) -> bool {
        match self {
            StateRange::M(range) => range.is_empty(),
            StateRange::D(range) => range.is_empty(),
            StateRange::I(range) => range.is_empty(),
        }
    }
}

/// The range within the query and reference where a continguous block of
/// matches occurs.
///
/// ## Validity
///
/// The two ranges should be the same length.
#[derive(Clone, Debug)]
pub struct MatchRange {
    /// The 0-based end-exclusive range of the match within the query.
    pub query_range: Range<usize>,
    /// The 0-based end-exclusive range of the match within the reference.
    pub ref_range:   Range<usize>,
}

impl MatchRange {
    /// Returns the length of the [`MatchRange`] in query and reference
    /// coordinates.
    pub fn len(&self) -> usize {
        // Validity: The ranges are the same length
        self.query_range.len()
    }

    /// Returns whether the [`MatchRange`] in empty in query and reference
    /// coordinates.
    pub fn is_empty(&self) -> bool {
        // Validity: The ranges are the same length
        self.query_range.is_empty()
    }
}

/// The range within the reference where a deletion occurs.
#[derive(Clone, Debug)]
pub struct DeletionRange {
    /// The 0-based end-exclusive range of the deletion within the reference.
    pub ref_range: Range<usize>,
}

impl DeletionRange {
    /// Returns the length of the [`DeletionRange`] in reference coordinates.
    ///
    /// The deletion has no length in query coordinates.
    pub fn len(&self) -> usize {
        self.ref_range.len()
    }

    /// Returns whether the [`DeletionRange`] in empty in reference coordinates.
    ///
    /// The deletion has no length in query coordinates.
    pub fn is_empty(&self) -> bool {
        // Validity: The ranges are the same length
        self.ref_range.is_empty()
    }
}

/// The range within the query where an insertion occurs, as well as the
/// corresponding index in the reference where it occurs.
#[derive(Clone, Debug)]
pub struct InsertionRange {
    /// The 0-based index of the insertion in the reference.
    pub ref_index:   InsertionIdx,
    /// The 0-based end-exclusive range of the insertion within the query.
    pub query_range: Range<usize>,
}

impl InsertionRange {
    /// Returns the length of the [`InsertionRange`] in query coordinates.
    ///
    /// The insertion has no length in reference coordinates.
    pub fn len(&self) -> usize {
        self.query_range.len()
    }

    /// Returns whether the [`InsertionRange`] in empty in query coordinates.
    ///
    /// The insertion has no length in reference coordinates.
    pub fn is_empty(&self) -> bool {
        self.query_range.is_empty()
    }
}

/// Alignment state ranges converted to coding sequence coordinates by
/// intersecting them with the exon ranges.
#[derive(Clone, Debug)]
pub enum CdsStateRange {
    M(CdsMatchRange),
    D(CdsDeletionRange),
    I(CdsInsertionRange),
}

impl CdsStateRange {
    pub fn match_range(&self) -> Option<&CdsMatchRange> {
        match self {
            CdsStateRange::M(range) => Some(range),
            _ => None,
        }
    }

    pub fn deletion_range(&self) -> Option<&CdsDeletionRange> {
        match self {
            CdsStateRange::D(range) => Some(range),
            _ => None,
        }
    }

    pub fn insertion_range(&self) -> Option<&CdsInsertionRange> {
        match self {
            CdsStateRange::I(range) => Some(range),
            _ => None,
        }
    }

    pub fn cds_len(&self) -> usize {
        match self {
            CdsStateRange::M(range) => range.cds_range.len(),
            CdsStateRange::D(range) => range.cds_range.len(),
            CdsStateRange::I(_) => 0,
        }
    }

    pub fn query_len(&self) -> usize {
        match self {
            CdsStateRange::M(range) => range.query_range.len(),
            CdsStateRange::D(_) => 0,
            CdsStateRange::I(range) => range.query_range.len(),
        }
    }
}

/// The range within the query and coding sequence where a contiguous block of
/// matches occurs.
///
/// ## Validity
///
/// The two ranges should be the same length.
#[derive(Clone, Debug)]
pub struct CdsMatchRange {
    /// The 0-based end-exclusive range of the match within the query.
    pub query_range: Range<usize>,
    /// The 0-based end-exclusive range of the match within the coding sequence.
    pub cds_range:   Range<usize>,
}

impl CdsMatchRange {
    /// Returns the lengths of the ranges.
    #[inline]
    pub(crate) fn len(&self) -> usize {
        // The ranges are the same length
        self.cds_range.len()
    }

    /// Returns whether the range is empty.
    #[inline]
    pub(crate) fn is_empty(&self) -> bool {
        // The ranges are the same length
        self.cds_range.is_empty()
    }

    /// Extends the start of both ranges by `amount`.
    ///
    /// This _decreases_ the start of the ranges.
    pub(crate) fn extend_start(&mut self, amount: usize) {
        // Validity: The ranges are extended by the same amount, so they remain
        // the same length
        self.cds_range.start -= amount;
        self.query_range.start -= amount;
    }

    /// Extends the end of both ranges by `amount`.
    ///
    /// This _increases_ the end of the ranges.
    pub(crate) fn extend_end(&mut self, amount: usize) {
        // Validity: The ranges are extended by the same amount, so they remain
        // the same length
        self.cds_range.end += amount;
        self.query_range.end += amount;
    }

    /// Cuts indices from the start of both ranges by `amount`.
    ///
    /// This _increases_ the start of the ranges.
    pub(crate) fn cut_start(&mut self, amount: usize) {
        // Validity: The ranges are cut by the same amount, so they remain the
        // same length
        self.cds_range.start += amount;
        self.query_range.start += amount;

        debug_assert!(self.cds_range.start <= self.cds_range.end);
        debug_assert!(self.query_range.start <= self.query_range.end);
    }

    /// Cuts indices from the end of both ranges by `amount`.
    ///
    /// This _decreases_ the end of the ranges.
    pub(crate) fn cut_end(&mut self, amount: usize) {
        // Validity: The ranges are cut by the same amount, so they remain the
        // same length
        self.cds_range.end -= amount;
        self.query_range.end -= amount;

        debug_assert!(self.cds_range.start <= self.cds_range.end);
        debug_assert!(self.query_range.start <= self.query_range.end);
    }
}

/// The range within the coding sequence where a deletion occurs.
#[derive(Clone, Debug)]
pub struct CdsDeletionRange {
    /// The range of the deletion within the coding sequence (non-empty,
    /// 0-based, end-exclusive).
    pub cds_range: Range<usize>,
}

impl CdsDeletionRange {
    /// Returns the length of the deletion in nucleotides.
    #[inline]
    pub(crate) fn len(&self) -> usize {
        self.cds_range.len()
    }
}

/// The range within the query where an insertion occurs, as well as the
/// corresponding index in the coding sequence where it occurs.
#[derive(Clone, Debug)]
pub struct CdsInsertionRange {
    /// The 0-based index of the insertion in the coding sequence.
    pub cds_index:   InsertionIdx,
    /// The 0-based end-exclusive range of the insertion within the query.
    pub query_range: Range<usize>,
}

impl CdsInsertionRange {
    /// Returns the length of the insertion.
    #[inline]
    pub(crate) fn len(&self) -> usize {
        self.query_range.len()
    }
}

impl StateRange {
    pub fn is_match(&self) -> bool {
        matches!(self, StateRange::M(_))
    }

    pub fn is_insert(&self) -> bool {
        matches!(self, StateRange::I(_))
    }

    pub fn is_delete(&self) -> bool {
        matches!(self, StateRange::D(_))
    }

    /// Returns the inclusive start index of the state in reference coordinates.
    ///
    /// For an insertion, this returns the index right of the insertion. When
    /// called on the first element of `Vec<StateRange>`, the right index will
    /// be the same as the start of a subsequent match or deletion.
    pub(crate) fn begin_ref_coord(&self) -> usize {
        match self {
            StateRange::M(state) => state.ref_range.start,
            StateRange::D(state) => state.ref_range.start,
            StateRange::I(state) => state.ref_index.right(),
        }
    }

    /// Returns the exclusive end index of the state in reference coordinates.
    ///
    /// For an insertion, this returns the index right of the insertion. When
    /// called on the last element of `Vec<StateRange>`, the right index will be
    /// the same as the exclusive end index of a preceding match or deletion.
    pub(crate) fn end_ref_coord(&self) -> usize {
        match self {
            StateRange::M(state) => state.ref_range.end,
            StateRange::D(state) => state.ref_range.end,
            // right is the exclusive end, whereas left would be inclusive
            StateRange::I(state) => state.ref_index.right(),
        }
    }

    /// Converts an [`Alignment`] to a sequence of [`StateRange`] for coordinate
    /// manipulation.
    ///
    /// All returned [`StateRange`] entries are guaranteed to be non-empty.
    ///
    /// ## Validity
    ///
    /// A valid `alignment` must be passed to ensure the output sequence of
    /// ranges will form an ordered partition of `query_range` and `ref_range`
    /// without overlap or gap. The alignment must only contain the operations
    /// `MIDS=X`.
    ///
    /// Furthermore, excluding soft clipping, the first/last states must be `M`,
    /// so that the output also has the first and last states as
    /// [`StateRange::M`].
    pub(crate) fn state_ranges_from_aligment<T>(alignment: &Alignment<T>) -> Vec<Self> {
        // This will be a slight overestimate due to the possible presence of
        // clipping.
        let mut states = Vec::with_capacity(alignment.states.len());

        // Initialize the current start of the next StateRange to the start of
        // the alignment within the query and the reference.
        let mut query_start = alignment.query_range.start;
        let mut ref_start = alignment.ref_range.start;

        // Validity: inc is nonzero based on guarantees of AlignmentStates, so
        // all StateRange entries will be non-empty
        for Ciglet { inc, op } in &alignment.states {
            match op {
                b'M' | b'=' | b'X' => {
                    // Validity: both ranges are length inc
                    states.push(Self::M(MatchRange {
                        query_range: query_start..query_start + inc,
                        ref_range:   ref_start..ref_start + inc,
                    }));
                    query_start += inc;
                    ref_start += inc;
                }
                b'I' => {
                    states.push(Self::I(InsertionRange {
                        // The insertion comes before ref_start
                        ref_index:   InsertionIdx::from_right_idx(ref_start),
                        query_range: query_start..query_start + inc,
                    }));
                    query_start += inc
                }
                b'D' => {
                    states.push(Self::D(DeletionRange {
                        ref_range: ref_start..ref_start + inc,
                    }));
                    ref_start += inc;
                }
                _ => {}
            }
        }

        states
    }
}

/// A coding sequence coordinate, which is either a [`Range`] for a match state
/// or an [`InsertionIdx`] where an insertion occurs.
#[derive(Clone, Eq, PartialEq, Hash, Debug)]
pub enum CdsCoord {
    /// A range of the coding sequence which is matched.
    M(Range<usize>),
    /// The index within the coding sequence where an insertion occurred.
    I(InsertionIdx),
}

/// An extension trait for basic 0-based range manipulation.
pub(crate) trait RangeExt: Sized {
    /// Adds a constant value to the range, shifting it right.
    #[must_use]
    fn add(&self, n: usize) -> Self;

    /// Subtracts a constant value from the range, shifting it left.
    #[must_use]
    fn sub(&self, n: usize) -> Self;

    /// Computes the shrinkage that was applied to go from `self` to `subset`.
    fn compute_shrinkage(&self, subset: &Self) -> RangeShrinkage;

    /// Shrinks the range by the given [`RangeShrinkage`] amount.
    ///
    /// The amount can be computed from [`compute_shrinkage`] or manually
    /// defined.
    ///
    /// [`compute_shrinkage`]: RangeExt::compute_shrinkage
    fn shrink(&self, amount: RangeShrinkage) -> Self;

    /// Checks whether the range overlaps with another.
    #[must_use]
    #[allow(dead_code)]
    fn overlaps(&self, other: &Self) -> bool;

    /// Checks whether the range contains an insertion strictly inside of it.
    #[must_use]
    fn contains_ins(&self, ins: InsertionIdx) -> bool;

    /// Intersects the two ranges, returning `Some` only if the result is
    /// non-empty.
    #[must_use]
    fn intersect(&self, other: &Self) -> Option<Self>;

    /// Intersects the end of `self` with the beginning of `next`, returning
    /// `Some` only if the result is non-empty.
    #[must_use]
    fn intersect_ordered(&self, next: &Self) -> Option<Self>;

    /// Checks whether the range is a superset of `other`.
    ///
    /// If the ranges are equal, this returns `true` (it is not strict).
    #[must_use]
    fn is_superset_of(&self, other: &Self) -> bool;

    /// Compares the positions of the ranges, returning `None` if the ranges
    /// overlap but are not equal.
    #[must_use]
    #[allow(dead_code)]
    fn strict_cmp(&self, other: &Self) -> Option<Ordering>;

    /// Compares the positions of the ranges, returning `None` if one range
    /// contains the other as a strict subset. Overlap is permitted.
    #[must_use]
    fn relaxed_cmp(&self, other: &Self) -> Option<Ordering>;

    /// Compare the position of the range with the insertion, returning
    /// [`Ordering::Less`] if the range comes before the insertion,
    /// [`Ordering::Greater`] if the range comes after, or [`None`] if the range
    /// strictly contains the insertion.
    fn cmp_ins(&self, ins: &InsertionIdx) -> Option<Ordering>;
}

impl RangeExt for Range<usize> {
    fn add(&self, n: usize) -> Self {
        self.start + n..self.end + n
    }

    fn sub(&self, n: usize) -> Self {
        self.start - n..self.end - n
    }

    fn compute_shrinkage(&self, subset: &Self) -> RangeShrinkage {
        debug_assert!(self.is_superset_of(subset), "The range must be a subset of `self`!");

        let shrink_left = subset.start - self.start;
        let shrink_right = self.end - subset.end;

        RangeShrinkage {
            shrink_left,
            shrink_right,
        }
    }

    fn shrink(&self, amount: RangeShrinkage) -> Self {
        self.start + amount.shrink_left..self.end - amount.shrink_right
    }

    fn overlaps(&self, other: &Self) -> bool {
        dbg_check_endpoints(self);
        dbg_check_endpoints(other);

        // Both ranges must end strictly after the other one starts in order for
        // overlap to occur
        self.end > other.start && other.end > self.start
    }

    fn contains_ins(&self, ins: InsertionIdx) -> bool {
        // Intuitively, the range occurs in [start, end) or [start, end-1]. The
        // insertion occurs between indices at left+0.5 or right-0.5. Overlap
        // requires start < right-0.5 < end-1 (without concern for equality,
        // since right-0.5 is not an integer). start < right-0.5 if and only if
        // start < right. right-0.5 < end-1 if and only if right < end-0.5, or
        // right < end.
        self.start < ins.right() && ins.right() < self.end
    }

    fn intersect(&self, other: &Self) -> Option<Self> {
        dbg_check_endpoints(self);
        dbg_check_endpoints(other);

        let start = self.start.max(other.start);
        let end = self.end.min(other.end);
        (start < end).then_some(start..end)
    }

    fn intersect_ordered(&self, next: &Self) -> Option<Self> {
        dbg_check_endpoints(self);
        dbg_check_endpoints(next);

        (next.start < self.end).then_some(next.start..self.end)
    }

    fn is_superset_of(&self, other: &Self) -> bool {
        dbg_check_endpoints(self);
        dbg_check_endpoints(other);

        self.start <= other.start && other.end <= self.end
    }

    fn strict_cmp(&self, other: &Self) -> Option<Ordering> {
        dbg_check_endpoints(self);
        dbg_check_endpoints(other);

        if self == other {
            Some(Ordering::Equal)
        } else if self.end <= other.start {
            Some(Ordering::Less)
        } else if self.start >= other.end {
            Some(Ordering::Greater)
        } else {
            None
        }
    }

    fn relaxed_cmp(&self, other: &Self) -> Option<Ordering> {
        dbg_check_endpoints(self);
        dbg_check_endpoints(other);

        match (self.start.cmp(&other.start), self.end.cmp(&other.end)) {
            // Handle equality first
            (Ordering::Equal, Ordering::Equal) => Some(Ordering::Equal),

            // At least one is strictly less per first case failing
            (Ordering::Less | Ordering::Equal, Ordering::Less | Ordering::Equal) => Some(Ordering::Less),

            // At least one is strictly greater per first case failing
            (Ordering::Greater | Ordering::Equal, Ordering::Greater | Ordering::Equal) => Some(Ordering::Greater),

            // One range contains the other
            (Ordering::Less, Ordering::Greater) | (Ordering::Greater, Ordering::Less) => None,
        }
    }

    fn cmp_ins(&self, ins: &InsertionIdx) -> Option<Ordering> {
        dbg_check_endpoints(self);

        if ins.right() <= self.start {
            // The index right of the insertion is the start of the range, or
            // less. The range is thus right of the insertion
            Some(Ordering::Greater)
        } else if self.end <= ins.right() {
            // The end-exclusive index is at most the index right of the
            // insertion. Equivalently, the end-inclusive index is at most the
            // index left of the insertion. The range is thus left of the
            // insertion.
            Some(Ordering::Less)
        } else {
            // This is precisely the case checked for in contains_ins
            None
        }
    }
}

pub(crate) struct RangeShrinkage {
    shrink_left:  usize,
    shrink_right: usize,
}

/// A wrapper around a coordinate-related type such that it implements 1-based,
/// end-inclusive [`Display`]. Semicolons are used as delimiters, and the same
/// `..` syntax is used for ranges rather than `..=`.
pub(crate) struct InclusiveFormatter<T>(T);

impl Display for InclusiveFormatter<&Range<usize>> {
    #[inline]
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, "{}..{}", self.0.start + 1, self.0.end)
    }
}

impl Display for InclusiveFormatter<&InsertionIdx> {
    #[inline]
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0.left_pos())
    }
}

impl Display for InclusiveFormatter<&CdsCoord> {
    #[inline]
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        match &self.0 {
            CdsCoord::M(range) => write!(f, "{}", range.display_inclusive()),
            CdsCoord::I(ins) => write!(f, "{}", ins.display_inclusive()),
        }
    }
}

impl<'a, T> Display for InclusiveFormatter<&'a [T]>
where
    InclusiveFormatter<&'a T>: Display,
{
    #[inline]
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        let Some((first, rest)) = self.0.split_first() else {
            return Ok(());
        };

        write!(f, "{}", InclusiveFormatter(first))?;

        rest.iter().try_for_each(|range| write!(f, ";{}", InclusiveFormatter(range)))
    }
}

impl<T> Display for InclusiveFormatter<&Nullable<T>>
where
    for<'a> InclusiveFormatter<&'a T>: Display,
    T: NullableValue,
{
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        match self.0.as_option() {
            Some(val) => write!(f, "{}", InclusiveFormatter(val)),
            None => Ok(()),
        }
    }
}

/// A trait for easily constructing [`InclusiveFormatter`].
pub(crate) trait InclusiveDisplay: Sized {
    /// Gets a 1-based end-inclusive display for a coordinate-related type.
    /// Semicolons are used as delimiters, and the same `..` syntax is used for
    /// ranges rather than `..=`.
    #[inline]
    #[must_use]
    fn display_inclusive(self) -> InclusiveFormatter<Self> {
        InclusiveFormatter(self)
    }
}

impl<T> InclusiveDisplay for T where for<'a> InclusiveFormatter<T>: Display {}

/// A trait implemented for types that can be parsed from a one-based string
/// representation.
///
/// This provides the opposite functionality as [`InclusiveDisplay`].
pub(crate) trait ParseOneBasedInclusive: Sized {
    /// Parses the one-based string representation of the coordinates.
    ///
    /// ## Errors
    ///
    /// See the implementations for possible parsing errors.
    fn parse_inclusive(coords: &str) -> std::io::Result<Self>;
}

impl ParseOneBasedInclusive for Range<usize> {
    /// Parses the one-based string representation of the coordinates.
    ///
    /// ## Errors
    ///
    /// The range must contain exactly one instance of `..`. The starting and
    /// ending coordinates must successfully parse to `usize`. The ending
    /// coordinate must be at least the starting coordinate, and the starting
    /// coordinate must be at least 1.
    fn parse_inclusive(coords: &str) -> std::io::Result<Self> {
        let coords = coords.trim();
        let mut range_parts = coords.split("..");
        let (Some(start), Some(end), None) = (range_parts.next(), range_parts.next(), range_parts.next()) else {
            return Err(std::io::Error::other(format!("Invalid coordinate range '{coords}'")));
        };
        parse_coordinate_range_from_parts(start, end)
    }
}

impl ParseOneBasedInclusive for InsertionIdx {
    /// Parses the one-based string representation of the insertion index
    /// coordinate.
    ///
    /// ## Errors
    ///
    /// The string must successfully parse as a non-zero `usize`.
    fn parse_inclusive(coords: &str) -> std::io::Result<Self> {
        let left_pos: usize = coords.parse().with_context(format!("Invalid position '{coords}'"))?;

        if left_pos == 0 {
            return Err(std::io::Error::other("Position must be at least 1"));
        }

        Ok(InsertionIdx::from_left_pos(left_pos))
    }
}

impl ParseOneBasedInclusive for CdsCoord {
    /// Parses the one-based string representation of the coordinates.
    ///
    /// ## Errors
    ///
    /// If the range does not contain `..`, it must successfully parse as a
    /// non-zero `usize`.
    ///
    /// Otherwise, it must contain exactly one instance of `..`, separating a
    /// starting and ending coordinate. The starting and ending coordinates must
    /// successfully parse to `usize`. The ending coordinate must be at least
    /// the starting coordinate, and the starting coordinate must be at least 1.
    fn parse_inclusive(coords: &str) -> std::io::Result<Self> {
        let coords = coords.trim();
        let mut range_parts = coords.split("..");
        let (Some(start), end, None) = (range_parts.next(), range_parts.next(), range_parts.next()) else {
            return Err(std::io::Error::other(format!("Invalid coordinates '{coords}'")));
        };

        match end {
            Some(end) => parse_coordinate_range_from_parts(start, end).map(CdsCoord::M),
            None => InsertionIdx::parse_inclusive(start).map(CdsCoord::I),
        }
    }
}

impl<T: ParseOneBasedInclusive> ParseOneBasedInclusive for Vec<T> {
    /// Parses the one-based string representation of the semicolon-delimited
    /// coordinates.
    ///
    /// If the string is empty or equal to [`HADOOP_NULL`], then the vector will
    /// be empty.
    ///
    /// ## Errors
    ///
    /// Each semicolon-delimited value must successfully parse. See the
    /// implementation of [`parse_inclusive`] for `T`.
    ///
    /// [`parse_inclusive`]: ParseOneBasedInclusive::parse_inclusive
    fn parse_inclusive(coords: &str) -> std::io::Result<Self> {
        if coords.is_empty() {
            return Ok(Vec::new());
        }

        coords.split(';').map(T::parse_inclusive).collect()
    }
}

/// An iterator for parsing semicolon-delimited one-based inclusive coordinates.
///
/// This is a lazy version of [`ParseOneBasedInclusive`] as implemented on
/// `Vec`. It is constructed with [`parse_coords_inclusive`].
///
/// ## Parameters
///
/// - `'a`: The lifetime of the string being parsed.
/// - `T`: The type to parse the coordinates into.
pub(crate) struct OneBasedInclusiveCoordsParser<'a, T> {
    /// The semicolon-delimited parts.
    parts:   std::str::Split<'a, char>,
    /// The type being parsed into.
    phantom: PhantomData<T>,
}

impl<T: ParseOneBasedInclusive> Iterator for OneBasedInclusiveCoordsParser<'_, T> {
    type Item = std::io::Result<T>;

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        self.parts.next().map(ParseOneBasedInclusive::parse_inclusive)
    }
}

/// Forms an iterator for lazily parsing the one-based inclusive ranges in a
/// semicolon-delimited string.
///
/// If the string is empty, then the iterator will be empty.
pub(crate) fn parse_coords_inclusive<T>(coords: &str) -> OneBasedInclusiveCoordsParser<'_, T> {
    let mut parts = coords.split(';');

    if coords.is_empty() {
        parts.next();
    }

    OneBasedInclusiveCoordsParser {
        parts,
        phantom: PhantomData,
    }
}

/// Given a 1-based starting coordinate and an ending coordinate as strings,
/// parses a non-empty 0-based range.
///
/// ## Errors
///
/// The `start` and `end` coordinates must successfully parse to `usize`. The
/// ending coordinate must be at least the starting coordinate, and the starting
/// coordinate must be at least 1.
fn parse_coordinate_range_from_parts(start: &str, end: &str) -> std::io::Result<Range<usize>> {
    // Parse 1-based inclusive range
    let start: usize = start.parse().with_context(format!("Invalid start coordinate '{start}'"))?;
    let end: usize = end.parse().with_context(format!("Invalid end coordinate '{end}'"))?;

    // Since we are using inclusive range, this also requires the range is
    // non-empty
    if end < start {
        return Err(std::io::Error::other(format!(
            "End coordinate must be >= start ({start}..{end})",
        )));
    }

    // Convert to 0-based half-open range (inclusive start, exclusive end)
    let Some(start) = start.checked_sub(1) else {
        return Err(std::io::Error::other(
            "Start coordinate must be at least 1 (for 1-based inclusive ranges)",
        ));
    };

    Ok(start..end)
}

/// When in debug mode, checks to make sure the end of the range is not less
/// than the start.
fn dbg_check_endpoints(range: &Range<usize>) {
    debug_assert!(
        range.start <= range.end,
        "The end of the range cannot be less than the start! Found {}..{}",
        range.start,
        range.end
    );
}
