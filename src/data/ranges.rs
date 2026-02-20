use crate::data::exons::ExonCoords;
use std::ops::Range;
use zoe::{alignment::Alignment, data::cigar::Ciglet};

/// Helper struct for protein annotationa and coordinate manipulation.
#[derive(Clone, Debug)]
pub(crate) enum StateRange {
    M(MatchRanges),
    D(DeletionRange),
    I(InsertionRange),
}

#[derive(Clone, Debug)]
pub(crate) struct MatchRanges {
    pub(crate) query_range: Range<usize>,
    pub(crate) ref_range:   Range<usize>,
}

#[derive(Clone, Debug)]
pub(crate) struct DeletionRange {
    pub(crate) ref_range: Range<usize>,
}

#[derive(Clone, Debug)]
pub(crate) struct InsertionRange {
    pub(crate) upstream_ref_index: usize,
    pub(crate) query_range:        Range<usize>,
}

impl MatchRanges {
    /// Only shift the `query_range` right (add) and not the reference.
    pub(crate) fn shift_query_right(&mut self, amount: usize) {
        self.query_range = self.query_range.start + amount..self.query_range.end + amount;
    }
}

impl InsertionRange {
    /// Only shift the `query_range` right (add) and not the CDS.
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

    /// Consolidate adjacent deletion states
    pub(crate) fn merge(&mut self, other: &Self) {
        debug_assert_eq!(self.cds_range.end, other.cds_range.start);

        self.cds_range.end = other.cds_range.end;
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

impl MatchRanges {
    fn intersect_exon(&self, exon: &ExonCoords) -> Option<CdsMatchRanges> {
        if self.ref_range.end <= exon.ref_range.start || exon.ref_range.end <= self.ref_range.start {
            None
        } else {
            let end_diff = self.ref_range.end.saturating_sub(exon.ref_range.end);
            let start_diff = exon.ref_range.start.saturating_sub(self.ref_range.start);
            let clipped_ref_start = self.ref_range.start + start_diff;
            let clipped_ref_end = self.ref_range.end - end_diff;
            Some(CdsMatchRanges {
                query_range: self.query_range.start + start_diff..self.query_range.end - end_diff,
                cds_range:   clipped_ref_start - exon.ref_to_cds_offset..clipped_ref_end - exon.ref_to_cds_offset,
            })
        }
    }
}

impl DeletionRange {
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

    pub(crate) fn state_ranges_from_aligment<T>(alignment: &Alignment<T>) -> Vec<Self> {
        let mut states = Vec::with_capacity(alignment.states.len());
        let mut query_start = alignment.query_range.start;
        let mut ref_start = alignment.ref_range.start;
        for Ciglet { inc, op } in &alignment.states {
            match op {
                b'M' | b'=' | b'X' => {
                    states.push(Self::M(MatchRanges {
                        query_range: query_start..query_start + inc,
                        ref_range:   ref_start..ref_start + inc,
                    }));
                    query_start += inc;
                    ref_start += inc;
                }
                b'I' => {
                    states.push(Self::I(InsertionRange {
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
                // Soft clipping is included in the ranges
                _ => {}
            }
        }

        states
    }
}
