//! Structs and enums for representing product and genome alignments,
//! specifically the ranges where matches, insertions, and deletions apply.

use crate::data::exons::ExonCoords;
use std::{
    cmp::Ordering,
    fmt::{self, Display, Formatter},
    ops::Range,
};
use zoe::{alignment::Alignment, data::cigar::Ciglet};

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

    /// The number of amino acids within the insertion's codon that are left of
    /// the insertion.
    ///
    /// This should be called only on [`InsertionIdx`] within a nucleotides
    /// sequence (not an amino acids sequence).
    ///
    /// 0 means that the insertion occurs in-frame. 1 means that the insertion
    /// appears after the first base of a codon. 2 means that the insertion
    /// appears after the second base of a codon.
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

/// The range within the reference where a deletion occurs.
#[derive(Clone, Debug)]
pub struct DeletionRange {
    /// The 0-based end-exclusive range of the deletion within the reference.
    pub ref_range: Range<usize>,
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

/// Alignment state ranges converted to coding sequence coordinates by
/// intersecting them with the exon ranges.
#[derive(Clone, Debug)]
pub enum CdsStateRange {
    M(CdsMatchRange),
    D(CdsDeletionRange),
    I(CdsInsertionRange),
}

impl CdsStateRange {
    /// Extracts a mutable reference to a [`CdsMatchRange`], or `None` is a
    /// different variant is present.
    pub fn match_range_mut(&mut self) -> Option<&mut CdsMatchRange> {
        match self {
            CdsStateRange::M(cds_match_ranges) => Some(cds_match_ranges),
            _ => None,
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
    /// Shifts the deletion in the coding sequence to the left by `amount`.
    ///
    /// This subtracts `amount` from the start and end of the range.
    pub(crate) fn shift_left(&mut self, amount: usize) {
        self.cds_range = self.cds_range.sub(amount);
    }

    /// Shifts the deletion in the coding sequence to the right by `amount`.
    ///
    /// This adds `amount` to the start and end of the range.
    pub(crate) fn shift_right(&mut self, amount: usize) {
        self.cds_range = self.cds_range.add(amount);
    }

    /// Returns the length of the deletion.
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
    // TODO: What is this used for? Why are both fields modified?
    /// Shift state left (subtract) by offset
    pub(crate) fn shift_left(&mut self, amount: usize) {
        self.cds_index.right -= amount;
        self.query_range = self.query_range.start - amount..self.query_range.end - amount;
    }

    // TODO: What is this used for? Why are both fields modified?
    /// Shift state right (add) by offset
    pub(crate) fn shift_right(&mut self, amount: usize) {
        self.cds_index.right += amount;
        self.query_range = self.query_range.start + amount..self.query_range.end + amount;
    }

    /// Returns the length of the insertion.
    pub(crate) fn len(&self) -> usize {
        self.query_range.len()
    }
}

impl MatchRange {
    fn intersect_exon(&self, exon: &ExonCoords) -> Option<CdsMatchRange> {
        // Intersect the two ranges in reference coordinates
        self.ref_range.intersect(&exon.ref_range).map(|intersect_ref_range| {
            let self_shrinkage = self.ref_range.compute_shrinkage(&intersect_ref_range);
            let intersect_query_range = self.query_range.shrink(self_shrinkage);

            // self.query_range and self.ref_range have the same length, so the
            // same should be true for the intersected versions
            debug_assert_eq!(intersect_query_range.len(), intersect_ref_range.len());

            let exon_shrinkage = exon.ref_range.compute_shrinkage(&intersect_ref_range);
            let intersect_cds_range = exon.cds_range.shrink(exon_shrinkage);

            // exon.cds_range and exon.ref_range have the same length, so the
            // same should be true for the intersected versions
            debug_assert_eq!(intersect_cds_range.len(), intersect_ref_range.len());

            // Validity: Per above, these are the same length, equal to the
            // intersect_ref_range.len()
            CdsMatchRange {
                query_range: intersect_query_range,
                cds_range:   intersect_cds_range,
            }
        })
    }
}

impl DeletionRange {
    fn intersect_exon(&self, exon: &ExonCoords) -> Option<CdsDeletionRange> {
        self.ref_range.intersect(&exon.ref_range).map(|intersect_ref_range| {
            let exon_shrinkage = exon.ref_range.compute_shrinkage(&intersect_ref_range);
            let intersect_cds_range = exon.cds_range.shrink(exon_shrinkage);

            CdsDeletionRange {
                cds_range: intersect_cds_range,
            }
        })
    }
}

impl InsertionRange {
    /// If the insertion range and exon strictly intersect (the insertion
    /// appears in the middle of the exon), then compute the
    /// [`CdsInsertionRange`] of the intersection.
    fn intersect_exon(&self, exon: &ExonCoords) -> Option<CdsInsertionRange> {
        // Intuitively, the range occurs in [start, end) or [start, end-1]. The
        // insertion occurs between indices at left+0.5 or right-0.5. Overlap
        // requires start < right-0.5 < end-1 (without concern for equality,
        // since right-0.5 is not an integer). start < right-0.5 if and only if
        // start < right. right-0.5 < end-1 if and only if right < end-0.5, or
        // right < end.

        (exon.ref_range.start < self.ref_index.right && self.ref_index.right < exon.ref_range.end).then(|| {
            // Intuitively, cds_range.start + (ref_index.right-ref_range.start)
            // where the second term is the offsetof the insertion within the
            // reference range. However, that offset may be positive or
            // negative, so we need to do additions before substractions to
            // prevent underflow
            let cds_index = InsertionIdx::from_right_idx(self.ref_index.right + exon.cds_range.start - exon.ref_range.start);

            CdsInsertionRange {
                cds_index,
                query_range: self.query_range.clone(),
            }
        })
    }
}

impl StateRange {
    /// Intersects a [`StateRange`] with an exon in reference coordinates,
    /// returning the resulting query coordinates and coding sequence
    /// coordinates.
    pub(crate) fn intersect_exon(&self, exon: &ExonCoords) -> Option<CdsStateRange> {
        match self {
            Self::M(m) => m.intersect_exon(exon).map(CdsStateRange::M),
            Self::D(d) => d.intersect_exon(exon).map(CdsStateRange::D),
            Self::I(i) => i.intersect_exon(exon).map(CdsStateRange::I),
        }
    }

    /// Converts an [`Alignment`] to a sequence of [`StateRange`] for coordinate
    /// manipulation.
    ///
    /// ## Validity
    ///
    /// A valid `alignment` must be passed to ensure the output sequence of
    /// ranges will form an ordered partition of `query_range` and `ref_range`
    /// without overlap or gap. In particular, this also requires the `N`
    /// operation to not be present in the `states`, which would otherwise cause
    /// a gap not covered in the reference.
    ///
    /// Furthermore, excluding soft clipping, the first/last `states` must be
    /// `M`, so that the output also has the first and last states as
    /// [`StateRange::M`].
    pub(crate) fn state_ranges_from_aligment<T>(alignment: &Alignment<T>) -> Vec<Self> {
        // This will be a slight overestimate due to the possible presence of
        // clipping.
        let mut states = Vec::with_capacity(alignment.states.len());

        // Initialize the current start of the next StateRange to the start of
        // the alignment within the query and the reference.
        let mut query_start = alignment.query_range.start;
        let mut ref_start = alignment.ref_range.start;

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
                b'N' => {
                    // TODO: It would be good not to support N
                    ref_start += inc;
                }
                // Soft clipping is included in the ranges, so no need to handle
                // S
                _ => {
                    // TODO: Add warning?
                }
            }
        }

        states
    }
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

    /// Intersects the two ranges, returning `Some` only if the result is
    /// non-empty.
    #[must_use]
    fn intersect(&self, other: &Self) -> Option<Self>;

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

    /// Returns a formatter displaying the 0-based range as 1-based, and
    /// end-inclusive instead of end-exclusive.
    ///
    /// See [`InclusiveRangeDisplay`] for more details.
    #[must_use]
    fn display_inclusive(&self) -> InclusiveRangeDisplay<'_>;
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

    fn intersect(&self, other: &Self) -> Option<Self> {
        dbg_check_endpoints(self);
        dbg_check_endpoints(other);

        let start = self.start.max(other.start);
        let end = self.end.min(other.end);
        (start < end).then_some(start..end)
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

    fn display_inclusive(&self) -> InclusiveRangeDisplay<'_> {
        InclusiveRangeDisplay { range: self }
    }
}

pub(crate) struct RangeShrinkage {
    shrink_left:  usize,
    shrink_right: usize,
}

/// A wrapper around a [`Range`] such that the end is displayed as 1-based and
/// end-inclusive.
///
/// Note that the same `..` syntax is used, rather than `..=`.
pub(crate) struct InclusiveRangeDisplay<'a> {
    /// The 0-based, end-exclusive range.
    range: &'a Range<usize>,
}

impl Display for InclusiveRangeDisplay<'_> {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, "{}..{}", self.range.start + 1, self.range.end)
    }
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
