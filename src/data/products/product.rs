use crate::{
    annotation::hashing::{nt_id, variant_hash},
    data::{
        products::{ComputedProduct, ProductSpec, incremental_products::ComputedIncrementalProducts},
        ranges::{CdsMatchRange, CdsStateRange, RangeExt},
    },
};
use std::{cmp::Ordering, ops::Range};
use zoe::prelude::*;

#[derive(Debug)]
pub(crate) struct Product<'a> {
    /// The ranges within the exons that the query covers. This is initially
    /// formed by intersecting the query ranges with the exon ranges, then is
    /// tweaked.
    ///
    /// This is guaranteed to contain ordered and non-overlapping ranges. It
    /// does not begin or end with [`CdsStateRange::I`]. Within the aligned
    /// portion of the coding sequence (the exons), the ranges will be adjacent
    /// (forming a partition). However, there may be exons at the beginning or
    /// end which are not aligned against (or partially aligned against).
    pub(crate) product_ranges:             Vec<CdsStateRange>,
    pub(crate) product_spec:               &'a ProductSpec,
    /// If this product's last exon ends at the stop extension position, this
    /// holds the query range of the stop extension nucleotides.
    pub(crate) stop_extension_query_range: Option<Range<usize>>,
}

impl<'a> Product<'a> {
    /// Returns true only if the requirement exists and the codon did not match.
    pub(crate) fn missing_required_start(&self, query: impl AsRef<[u8]>) -> bool {
        if let Some(required) = self.product_spec.exons.required_start
            && let Some(CdsStateRange::M(m)) = self.product_ranges.first()
            && m.cds_range.start == 0
            && m.cds_range.end >= 3
            && let Some(mut first) = query
                .as_ref()
                .get(m.query_range.start..)
                .and_then(|s| s.first_chunk().copied())
        {
            first.make_ascii_uppercase();
            first != required
        } else {
            false
        }
    }

    pub(crate) fn leading_cds_gap_len(&self) -> usize {
        self.product_ranges
            .iter()
            .find_map(|s| match s {
                CdsStateRange::M(m) => Some(m.cds_range.start),
                CdsStateRange::D(d) => Some(d.cds_range.start),
                _ => None,
            })
            .unwrap_or(0)
    }

    /// Because the exon intersection algorithm can split deletions spanning
    /// multiple exons, this method merges adjacent deletion CDS for downstream
    /// correctness.
    pub(crate) fn condense_deletions(&mut self) {
        // The order of the closure arguments in dedup_by gets reversed.
        self.product_ranges.dedup_by(|range2, range1| {
            if let CdsStateRange::D(del1) = range1
                && let CdsStateRange::D(del2) = range2
            {
                // Ensure the deletions are adjacent in the coding sequence.
                // Validity: product_ranges forms an ordered partition over the
                // portion of the coding sequence which is aligned against.
                // Hence, adjacent deletions in product_ranges are guaranteed to
                // also be adjacent in the coding sequence.
                debug_assert_eq!(del1.cds_range.end, del2.cds_range.start);
                // Extend del1 (the kept element) to encompass del2
                del1.cds_range.end = del2.cds_range.end;
                // Return true to signal that deduplication is needed
                true
            } else {
                false
            }
        });
    }

    /// Fix indel frames by repositioning insertions and deletions to in-frame
    /// boundaries.
    ///
    /// This corrects in-frame indels (length divisible by 3) that occur at
    /// out-of-frame positions by shifting them left or right based on codon
    /// usage statistics. The algorithm mirrors the logic in
    /// `codonCorrectStats.pl`, albeit for coordinate math instead of strings.
    ///
    /// ## Validity
    ///
    /// `query_seq` should contain unaligned, uppercase IUPAC bases.
    pub(crate) fn fix_frames(&mut self, query: impl AsRef<[u8]>) {
        let len = self.product_ranges.len();
        let query_seq = query.as_ref();

        // We iterate by index because we need to access neighbors
        let mut i = 0;
        while i < len {
            match &self.product_ranges[i] {
                CdsStateRange::I(ins) => {
                    // Upstream index becomes upstream 1-based position
                    let frame = ins.frame();

                    // Only correct in-frame insertions at out-of-frame positions
                    if frame != 0 && ins.len() % 3 == 0 {
                        // Validity: this function requires the same
                        // requirements on query_seq
                        self.fix_insertion_frame(i, query_seq);
                    }
                }
                CdsStateRange::D(del) => {
                    // The inclusive 0-based index is the upstream 1-based position
                    let frame = del.cds_range.start % 3;

                    // Only correct in-frame deletions at out-of-frame positions
                    if frame != 0 && del.len() % 3 == 0 {
                        self.fix_deletion_frame(i, query_seq);
                    }
                }
                CdsStateRange::M(_) => {}
            }
            i += 1;
        }
    }

    // TODO: Can fixing frame mess up order of partition in case of flanking
    // deletion?

    /// Fixes an insertion at `idx` by shifting it to an in-frame position.
    ///
    /// Uses the A1/A2 shift logic from `codonCorrectStats.pl`:
    /// - A1 (frame 1): Compare codons formed by:
    ///   - 1 preceding base + first 2 insert bases
    ///   - last 1 insert base + 2 following
    ///+
    /// - A2 (frame 2): Compare codons formed by:
    ///   - 2 preceding bases + first 1 insert base
    ///   - last 2 insert bases + 1 following
    ///
    /// ## Validity
    ///
    /// `query_seq` should contain unaligned, uppercase IUPAC bases.
    fn fix_insertion_frame(&mut self, idx: usize, query_seq: &[u8]) -> Option<()> {
        let Product {
            product_ranges,
            product_spec,
            ..
        } = self;

        // Insertions cannot be split over exons, and in practice many scoring
        // schemes tend to avoid adjacent insertions and deletions, so in
        // practice left_match and right_match occur at idx-1 and idx+1
        // respectively. But in case either is a deletion, we perform a search.
        let Some((left_match, CdsStateRange::I(ins), right_match)) = Self::partition_states(product_ranges, idx) else {
            return None;
        };

        let frame = ins.frame();
        let codon_index = ins.codon_index();
        let insert_len = ins.len();

        // for weight lookup
        let codon_position = codon_index + 1;
        let insert_seq = &query_seq[ins.query_range.clone()];

        if frame == 2 {
            // A2: insertion after 2nd codon base
            //
            // A2/L2: shift insert left 2
            //   New codon: last 2 bases of insert + cp3
            //
            // A2/R1: shift insert right 1
            //   New codon: cp1 + cp2 + first 1 base of insert
            let cp1 = *query_seq.get(left_match.query_range.end - 2)?;
            let cp2 = *query_seq.get(left_match.query_range.end - 1)?;
            let cp3 = *query_seq.get(right_match.query_range.start)?;

            let a2l2 = [insert_seq[insert_len - 2], insert_seq[insert_len - 1], cp3];
            let a2r1 = [cp1, cp2, insert_seq[0]];

            // Validity: The codons are uppercase because they are derived from
            // bases in query_seq
            if product_spec.codon_left_ge_right(a2r1, a2l2, codon_position as u32) {
                // Insertion shifts right 1 for frame 2 split codon
                left_match.extend_end(1);
                ins.shift_right(1);
                right_match.cut_start(1);
            } else {
                // Insertion shifts left 2 for frame 2 split codon
                left_match.cut_end(2);
                ins.shift_left(2);
                right_match.extend_start(2);
            }
        } else if frame == 1 {
            // A1 insertion: insertion after 1st base of codon
            //
            // A1/L1: shift insert left 1
            //   New codon: last 1 base of insert + cp2 + cp3
            //
            // A1/R2: shift insert right 2
            //   New codon: cp1 + first 2 bases of insert

            let cp1 = *query_seq.get(left_match.query_range.end - 1)?;
            let cp2 = *query_seq.get(right_match.query_range.start)?;
            let cp3 = *query_seq.get(right_match.query_range.start + 1)?;

            let a1l1 = [insert_seq[insert_len - 1], cp2, cp3];
            let a1r2 = [cp1, insert_seq[0], insert_seq[1]];

            // Validity: The codons are uppercase because they are derived from
            // bases in query_seq
            if product_spec.codon_left_ge_right(a1l1, a1r2, codon_position as u32) {
                // Insertion shifts left 1 for frame 1 split codon
                left_match.cut_end(1);
                ins.shift_left(1);
                right_match.extend_start(1);
            } else {
                // Insertion shifts right 2 for frame 1 split codon
                left_match.extend_end(2);
                ins.shift_right(2);
                right_match.cut_start(2);
            }
        }

        Some(())
    }

    /// Fix a deletion at index `i` by shifting it to an in-frame position.
    fn fix_deletion_frame(&mut self, idx: usize, query: &[u8]) -> Option<()> {
        let Product {
            product_ranges,
            product_spec,
            ..
        } = self;

        let Some((left_match, CdsStateRange::D(del), right_match)) = Self::partition_states(product_ranges, idx) else {
            return None;
        };

        // Codon positions (1-based) at the boundaries of the deletion for table
        // checking
        let pos_left = (del.cds_range.start / 3) + 1;
        let pos_right = ((del.cds_range.end - 1) / 3) + 1;

        let frame = del.cds_range.start % 3;

        // Get the pivot codon bases that cross the deletion boundary
        if frame == 1 {
            // Frame 1 pivot codon: is last 1 base before gap + first 2 bases
            // after gap
            let cp1 = query.get(left_match.query_range.end - 1)?.to_ascii_uppercase();
            let cp2 = query.get(right_match.query_range.start)?.to_ascii_uppercase();
            let cp3 = query.get(right_match.query_range.start + 1)?.to_ascii_uppercase();

            let pivot = [cp1, cp2, cp3];

            // TODO: Likely want to compare with Ordering::is_gt instead, so
            // that ties resolve as right shift.

            // By default, we shift the deletion left for frame 1 (causing 1
            // match to shift right, rather than 2). Only if there is evidence
            // for shifting the deletion to right do we do that. pos_left being
            // better means shifting matches left is better, which means
            // shifting deletion right is better. Validity: The codon is
            // uppercase because it is derived from bases in query
            let shift_del_right = product_spec
                .compare_codon_positions(pos_left as u32, pos_right as u32, pivot)
                .is_some_and(Ordering::is_ge);

            if shift_del_right {
                // Shift the deletion right by 2, causing 2 matches to move left
                left_match.extend_end(2);
                del.shift_right(2);
                right_match.cut_start(2);
            } else {
                // Shift the deletion left by 1, causing 1 match to move right
                left_match.cut_end(1);
                del.shift_left(1);
                right_match.extend_start(1);
            }
        } else if frame == 2 {
            // Frame2 pivot codon is last 2 bases before gap + first 1 base
            // after gap.
            let cp1 = query.get(left_match.query_range.end - 2)?.to_ascii_uppercase();
            let cp2 = query.get(left_match.query_range.end - 1)?.to_ascii_uppercase();
            let cp3 = query.get(right_match.query_range.start)?.to_ascii_uppercase();

            let pivot = [cp1, cp2, cp3];

            // By default, we shift the deletion right for frame 2 (causing 1
            // match to shift left, rather than 2). Only if there is evidence
            // for a left shift do we do that. pos_right being better means
            // shifting matches right is better, which means shifting deletion
            // left is better. Validity: The codon is uppercase because it is
            // derived from bases in query
            let shift_del_left = product_spec
                .compare_codon_positions(pos_left as u32, pos_right as u32, pivot)
                .is_some_and(Ordering::is_lt);

            if shift_del_left {
                // Shift the deletion left by 2, causing 2 matches to move right
                left_match.cut_end(2);
                del.shift_left(2);
                right_match.extend_start(2);
            } else {
                // Shift the deletion right by 1, causing 1 match to move right
                left_match.extend_end(1);
                del.shift_right(1);
                right_match.cut_start(1);
            }
        }

        Some(())
    }

    /// Get the match ranges to modify and the indel
    fn partition_states(
        product_ranges: &mut [CdsStateRange], idx: usize,
    ) -> Option<(&mut CdsMatchRange, &mut CdsStateRange, &mut CdsMatchRange)> {
        let (left, mid, right) = product_ranges.split_around_mut(idx)?;

        let l = left.iter_mut().rev().filter_map(CdsStateRange::match_range_mut).next()?;
        let r = right.iter_mut().filter_map(CdsStateRange::match_range_mut).next()?;

        Some((l, mid, r))
    }

    pub(crate) fn add_query_coords(&mut self, offset: usize) {
        if offset == 0 {
            return;
        }

        for state in &mut self.product_ranges {
            match state {
                CdsStateRange::M(m) => m.query_range = m.query_range.add(offset),
                CdsStateRange::I(i) => i.query_range = i.query_range.add(offset),
                _ => {}
            }
        }
    }
}

/// An extension trait for slices, providing functionality specific to
/// DAIS-ribosome.
pub(crate) trait SliceExt<T> {
    /// Extracts a mutable index from a slice, along with the slice before it
    /// and the slice after it.
    fn split_around_mut(&mut self, index: usize) -> Option<(&mut [T], &mut T, &mut [T])>;
}

impl<T> SliceExt<T> for [T] {
    fn split_around_mut(&mut self, index: usize) -> Option<(&mut [T], &mut T, &mut [T])> {
        let (left, rest) = self.split_at_mut_checked(index)?;
        let (mid, right) = rest.split_first_mut()?;
        Some((left, mid, right))
    }
}

impl<'a> Product<'a> {
    /// Computes the output data for this product, materializing all ranges into
    /// sequences using `query`.
    ///
    /// ## Validity
    ///
    /// The `query` should contain unaligned, uppercase IUPAC bases.
    pub fn materialize(&self, query: &Nucleotides) -> ComputedProduct<'a> {
        // Compute all the fields that rely on incremental updates until the
        // first stop codon
        let incremental = ComputedIncrementalProducts::new(query, self);

        let ComputedIncrementalProducts {
            cds_aln,
            aa_aln,
            cds_seq,
            has_insertion,
            has_shift_indel,
            query_coords,
            cds_coords,
            insertions,
            deletions,
            trailing_cds_unaligned,
        } = incremental;

        // Form aa_seq by splicing insertions into aa_aln (and removing
        // deletions)
        let aa_seq = {
            let mut out = AminoAcids::new();

            let mut aa_aln_without_deletions = aa_aln.iter().filter(|&&b| b != b'-' && b != b'.').copied();

            // The number of amino acids consumed from aa_aln so far
            let mut num_consumed = 0;
            for insertion in &insertions {
                // 1-based index after which insertion occurs is equivalently
                // the count of the number of amino acids before the insertion.
                let num_to_consume = insertion.upstream_aa - num_consumed;

                out.extend(aa_aln_without_deletions.by_ref().take(num_to_consume));
                out.extend_from_slice(&insertion.inserted_residues);

                num_consumed += num_to_consume;
            }

            // Consume the rest of the amino acids after the last insertion
            out.extend(aa_aln_without_deletions);

            out
        };

        // Get hashes
        let cds_id = nt_id(&cds_seq);
        let variant_hash = variant_hash(&aa_seq);

        ComputedProduct {
            product_name: &self.product_spec.name,
            cds_seq,
            cds_aln,
            cds_id,
            aa_seq,
            aa_aln,
            variant_hash,
            has_insertion,
            has_shift_indel,
            query_coords,
            cds_coords,
            insertions,
            deletions,
            trailing_cds_unaligned,
        }
    }
}
