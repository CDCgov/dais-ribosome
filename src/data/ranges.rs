use crate::data::exons::ExonCoords;
use std::ops::Range;
use zoe::{alignment::Alignment, data::cigar::Ciglet};

/// The range/index within the query and reference where a contiguous block of
/// matches, deletions, or insertions occur.
///
/// This is a helper type for protein annotation and coordinate manipulation.
#[derive(Clone, Debug)]
pub(crate) enum StateRange {
    M(MatchRange),
    D(DeletionRange),
    I(InsertionRange),
}

impl StateRange {
    /// Shifts any query indices to the right (addition) without altering any
    /// reference indices.
    pub(crate) fn shift_query_right(&mut self, amount: usize) {
        match self {
            StateRange::M(m) => m.shift_query_right(amount),
            StateRange::I(ins) => ins.shift_query_right(amount),
            StateRange::D(_) => {}
        }
    }
}

/// The range within the query and reference where a continguous block of
/// matches occurs.
///
/// The two ranges will be the same length.
#[derive(Clone, Debug)]
pub(crate) struct MatchRange {
    /// The 0-based end-exclusive range of the match within the query.
    pub(crate) query_range: Range<usize>,
    /// The 0-based end-exclusive range of the match within the reference.
    pub(crate) ref_range:   Range<usize>,
}

/// The range within the reference where a deletion occurs.
#[derive(Clone, Debug)]
pub(crate) struct DeletionRange {
    /// The 0-based end-exclusive range of the deletion within the reference.
    pub(crate) ref_range: Range<usize>,
}

/// The range within the query where an insertion occurs, as well as the
/// corresponding position in the reference TODO (at which, before which??)
#[derive(Clone, Debug)]
pub(crate) struct InsertionRange {
    // TODO: could rename this to make it clear
    /// The in the reference *after* which the insertion occurs.
    pub(crate) upstream_ref_index: usize,
    pub(crate) query_range:        Range<usize>,
}

impl MatchRange {
    /// Shifts the `query_range` to the right (addition) without altering the
    /// range in the reference.
    pub(crate) fn shift_query_right(&mut self, amount: usize) {
        self.query_range = self.query_range.start + amount..self.query_range.end + amount;
    }
}

impl InsertionRange {
    /// Shifts the `query_range` to the right (addition) without altering the
    /// index in the reference.
    pub(crate) fn shift_query_right(&mut self, amount: usize) {
        self.query_range = self.query_range.start + amount..self.query_range.end + amount;
    }
}

/// Alignment state ranges converted to CDS coordinates after exon intersection.
#[derive(Clone, Debug)]
pub(crate) enum CdsStateRange {
    M(CdsMatchRanges),
    D(CdsDeletionRange),
    I(CdsInsertionRange),
}

impl CdsStateRange {
    /// Extracts a mutable reference to a [`CdsMatchRanges`], or `None` is a
    /// different variant is present.
    pub fn match_range_mut(&mut self) -> Option<&mut CdsMatchRanges> {
        match self {
            CdsStateRange::M(cds_match_ranges) => Some(cds_match_ranges),
            _ => None,
        }
    }
}

#[derive(Clone, Debug)]
pub(crate) struct CdsMatchRanges {
    pub(crate) query_range: Range<usize>,
    pub(crate) cds_range:   Range<usize>,
}

impl CdsMatchRanges {
    /// Extend the end of both ranges by `amount`.
    pub(crate) fn extend_end(&mut self, amount: usize) {
        self.cds_range.end += amount;
        self.query_range.end += amount;
    }

    /// Shrink from the end of both ranges by `amount`.
    pub(crate) fn shrink_end(&mut self, amount: usize) {
        self.cds_range.end -= amount;
        self.query_range.end -= amount;
    }

    /// Extend the start of both ranges earlier by `amount`.
    pub(crate) fn extend_start(&mut self, amount: usize) {
        self.cds_range.start -= amount;
        self.query_range.start -= amount;
    }

    /// Shrink from the start of both ranges by `amount` (move start later).
    pub(crate) fn shrink_start(&mut self, amount: usize) {
        debug_assert!(amount <= self.cds_range.len());
        debug_assert!(amount <= self.query_range.len());

        self.cds_range.start += amount;
        self.query_range.start += amount;
    }

    /// Shift ONLY the `query_range` (add) right by the `offset`
    pub(crate) fn shift_query_right(&mut self, offset: usize) {
        self.query_range = self.query_range.start + offset..self.query_range.end + offset;
    }
}

#[derive(Clone, Debug)]
pub(crate) struct CdsDeletionRange {
    pub(crate) cds_range: Range<usize>,
}

impl CdsDeletionRange {
    /// Shift deletion CDS left (subtract) by amount
    pub(crate) fn shift_left(&mut self, amount: usize) {
        self.cds_range = self.cds_range.start - amount..self.cds_range.end - amount;
    }

    /// Shift deletion CDS right (add) by amount
    pub(crate) fn shift_right(&mut self, amount: usize) {
        self.cds_range = self.cds_range.start + amount..self.cds_range.end + amount;
    }

    pub(crate) fn len(&self) -> usize {
        self.cds_range.len()
    }
}

#[derive(Clone, Debug)]
pub(crate) struct CdsInsertionRange {
    pub(crate) upstream_cds_index: usize,
    pub(crate) query_range:        Range<usize>,
}

impl CdsInsertionRange {
    /// Shift state left (subtract) by offset
    pub(crate) fn shift_left(&mut self, amount: usize) {
        self.upstream_cds_index -= amount;
        self.query_range = self.query_range.start - amount..self.query_range.end - amount;
    }

    /// Only shift the `query_range` right (add) and not the CDS.
    pub(crate) fn shift_query_right(&mut self, amount: usize) {
        self.query_range = self.query_range.start + amount..self.query_range.end + amount;
    }

    /// Shift state right (add) by offset
    pub(crate) fn shift_right(&mut self, amount: usize) {
        self.upstream_cds_index += amount;
        self.query_range = self.query_range.start + amount..self.query_range.end + amount;
    }

    pub(crate) fn len(&self) -> usize {
        self.query_range.len()
    }
}

impl MatchRange {
    fn intersect_exon(&self, exon: &ExonCoords) -> Option<CdsMatchRanges> {
        if self.ref_range.end <= exon.ref_range.start || exon.ref_range.end <= self.ref_range.start {
            None
        } else {
            // TODO: The saturating sub followed by addition/subtraction is a
            // convoluted way of taking the minimum/maximum

            // The number of bases that the match range extends past the end of
            // the exon on the right
            let end_diff = self.ref_range.end.saturating_sub(exon.ref_range.end);
            // The number bases that the match range extends past the end of the
            // exon on the left
            let start_diff = exon.ref_range.start.saturating_sub(self.ref_range.start);
            // The start of the intersected range, in the reference coordinates
            let clipped_ref_start = self.ref_range.start + start_diff;
            // The end of the intersected range, in the reference coordinates
            let clipped_ref_end = self.ref_range.end - end_diff;

            Some(CdsMatchRanges {
                query_range: self.query_range.start + start_diff..self.query_range.end - end_diff,
                cds_range:   clipped_ref_start - exon.ref_to_cds_offset..clipped_ref_end - exon.ref_to_cds_offset,
            })
        }
    }
}

impl DeletionRange {
    // TODO: This is identical to above
    fn intersect_exon(&self, exon: &ExonCoords) -> Option<CdsDeletionRange> {
        if self.ref_range.end <= exon.ref_range.start || exon.ref_range.end <= self.ref_range.start {
            None
        } else {
            let end_diff = self.ref_range.end.saturating_sub(exon.ref_range.end);
            let start_diff = exon.ref_range.start.saturating_sub(self.ref_range.start);
            let clipped_ref_start = self.ref_range.start + start_diff;
            let clipped_ref_end = self.ref_range.end - end_diff;
            Some(CdsDeletionRange {
                cds_range: clipped_ref_start - exon.ref_to_cds_offset..clipped_ref_end - exon.ref_to_cds_offset,
            })
        }
    }
}

impl InsertionRange {
    fn intersect_exon(&self, exon: &ExonCoords) -> Option<CdsInsertionRange> {
        if self.upstream_ref_index < exon.ref_range.start || self.upstream_ref_index >= exon.ref_range.end - 1 {
            None
        } else {
            Some(CdsInsertionRange {
                upstream_cds_index: self.upstream_ref_index - exon.ref_to_cds_offset,
                query_range:        self.query_range.clone(),
            })
        }
    }
}

impl StateRange {
    /// Intersect with an exon and convert to CDS coordinates.
    pub(crate) fn intersect_exon(&self, exon: &ExonCoords) -> Option<CdsStateRange> {
        match self {
            Self::M(m) => m.intersect_exon(exon).map(CdsStateRange::M),
            Self::D(d) => d.intersect_exon(exon).map(CdsStateRange::D),
            Self::I(i) => i.intersect_exon(exon).map(CdsStateRange::I),
        }
    }

    /// Converts an [`Alignment`] to a sequence of [`StateRange`] for coordinate
    /// manipulation.
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
                    states.push(Self::M(MatchRange {
                        query_range: query_start..query_start + inc,
                        ref_range:   ref_start..ref_start + inc,
                    }));
                    query_start += inc;
                    ref_start += inc;
                }
                b'I' => {
                    states.push(Self::I(InsertionRange {
                        // TODO: THIS WILL PANIC IF ALIGNMENT STARTED AT 0????
                        upstream_ref_index: ref_start - 1,
                        query_range:        query_start..query_start + inc,
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
                    ref_start += inc;
                }
                // Soft clipping is included in the ranges, so no need to handle
                // S
                _ => {}
            }
        }

        states
    }
}
