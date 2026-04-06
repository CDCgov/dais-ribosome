//! The core process of DAIS-ribosome for creating a range-based
//! [`RibosomeOutput`] annotation result from an input [`QueryRecord`].

use crate::{
    AlignmentStatesExt,
    config::annotation_module::{AnnotationModule, ReferenceGroup},
    data::{
        QueryRecord,
        ranges::{CdsMatchRange, CdsStateRange, InsertionIdx, InsertionRange, RangeExt, StateRange},
    },
    error::RibosomeError,
    outputs::{GenomeAndProductStates, Product, RibosomeOutput},
};
use std::{cmp::Ordering, ops::Range};
use zoe::{
    alignment::Alignment,
    data::{SanitizeBase, types::nucleotides::CodonExtension},
    prelude::*,
};

impl<'a> AnnotationModule<'a> {
    /// Processes a single query, returning [`RibosomeOutput`] containing all
    /// the genome alignments against the relevant references, as well as the
    /// protein products formed for each reference.  
    pub fn process(&self, query: QueryRecord) -> Result<RibosomeOutput<'_>, RibosomeError> {
        // Get the corresponding reference information for the compound type of
        // the query
        let Some(reference_data) = self.ctype_map.get(&query.ctype) else {
            return Err(RibosomeError::UnimplementedCtype(query.ctype.into()));
        };

        let mut states = Vec::with_capacity(reference_data.len());

        for ref_id_data in reference_data.iter() {
            let (query_ori_offset, chewed_query) = self.rule_chew_to_start(&query, ref_id_data);

            // TODO: Do we ever do revcomp alignment?

            // Get the alignment to the best reference
            let Some(mut genome_aln) = ref_id_data.best_alignment(&chewed_query) else {
                return Err(format!("Query '{}' could not be aligned to any reference", query.id).into());
            };

            genome_aln.query_range = genome_aln.query_range.add(query_ori_offset);
            genome_aln.query_len += query_ori_offset;
            genome_aln.states.prepend_soft_clip(query_ori_offset);

            // Extend the left and right side of the alignments
            self.rule_repairable_ends(&mut genome_aln);

            // Compute the stop extension
            let stop_extension = self.rule_stop_extension(&query, &genome_aln);

            // Validity: requirements met based on best_alignment guarantees
            let genome_aln_states = StateRange::state_ranges_from_aligment(&genome_aln);

            let mut products = Vec::with_capacity(ref_id_data.proteins.len());

            for product in &ref_id_data.proteins {
                // Validity: requirements met based on
                // state_ranges_from_aligment guarantees
                let mut product_ranges = product.make_product_ranges(&genome_aln_states);

                // Validity: the same `query_seq` is passed as was used to form
                // `genome_aln_states`
                if product_ranges.missing_required_start(&query) {
                    continue;
                }

                product_ranges.condense_deletions();
                product_ranges.fix_frames(&query);

                products.push(product_ranges);
            }

            // Push stop extension into every product whose last exon ends at
            // the extension's reference position
            if let Some(ext) = &stop_extension {
                for product in &mut products {
                    // Check whether the last exon ends at the same place the
                    // stop extension "ends" (the index before which it occurs).
                    if let Some(last_exon) = product.product_spec.exons.coords.last()
                        && last_exon.ref_range.end == ext.ref_index.right()
                    {
                        product.stop_extension_query_range = Some(ext.query_range.clone());
                    }
                }
            }

            let ref_len = ref_id_data.length;

            let leading_ref_unaligned = genome_aln_states
                .iter()
                .find_map(|s| match s {
                    StateRange::M(m) => Some(m.ref_range.start),
                    StateRange::D(d) => Some(d.ref_range.start),
                    _ => None,
                })
                .unwrap_or(0);

            let trailing_ref_unaligned = ref_len
                - genome_aln_states
                    .iter()
                    .rev()
                    .find_map(|s| match s {
                        StateRange::M(m) => Some(m.ref_range.end),
                        StateRange::D(d) => Some(d.ref_range.end),
                        _ => None,
                    })
                    .unwrap_or(0);

            states.push(GenomeAndProductStates {
                reference_id: &ref_id_data.reference_id,
                ref_len,
                genome_aln_states,
                leading_ref_unaligned,
                trailing_ref_unaligned,
                products,
                stop_extension_query_range: stop_extension.map(|ins| ins.query_range),
            });
        }

        Ok(RibosomeOutput {
            query,
            states,
            formatting: &self.data.formatting,
        })
    }

    /// If the alignment reaches the end of the reference but does not end in a
    /// stop codon as expected, then attempts to represent any unaligned bases
    /// at the tail of the query up until the first stop codon as an insertion.
    ///
    /// This insertion is called the stop extension. The stop codon that is
    /// searched for must be in-frame.
    fn rule_stop_extension(&self, query: &QueryRecord, genome_aln: &Alignment<u32>) -> Option<InsertionRange> {
        // Note: This contains uppercase IUPAC, possibly with either U or T
        let query_seq = &query.nucleotides;

        if self.data.rules.list_contig_stop_extension
            && genome_aln.uanligned_ref_tail() == 0
            && genome_aln.unaligned_query_tail() >= 3
            && let Some(last_aligned_codon) = query_seq.slice(genome_aln.aln_query_range()).get_tail_codon()
            // Do not extend past known stop codon
            && !last_aligned_codon.is_std_stop_codon()
            && let Some(stop_codon_index) = query_seq.slice(genome_aln.query_range.end..).find_next_aa_in_frame(b'*')
        {
            // The exclusive end of the last alignment range is the inclusive
            // start of the insertion
            let start_index = genome_aln.query_range.end;

            // Convert inclusive start of codon to exclusive end of codon
            let end_index = genome_aln.query_range.end + stop_codon_index + 3;

            Some(InsertionRange {
                // The insertion is before the end (ref_range.end)
                ref_index:   InsertionIdx::from_right_idx(genome_aln.ref_range.end),
                query_range: start_index..end_index,
            })
        } else {
            None
        }
    }

    // TODO: What if we have:
    // Query      ----ATG----------
    // Reference      --------------
    // Then this rule won't apply

    /// If the `chew_to_start` rule is enabled, then handle the case when the
    /// query is longer than the reference.
    ///
    /// If a start codon can be found, such that excluding anything before the
    /// start codon results in a query that is at least the length of the
    /// reference, then slice the query to contain the start codon and
    /// everything after it. The first value returned is the starting position
    /// of the returned slice.
    ///
    /// If any of these conditions fail to hold, no shrinking occurs, and the
    /// starting position of the returned slice is 0.
    fn rule_chew_to_start<'b>(
        &self, query: &'b QueryRecord, ref_id_data: &ReferenceGroup<'_>,
    ) -> (usize, NucleotidesView<'b>) {
        if self.data.rules.chew_to_start
            && query.nucleotides.len() > ref_id_data.length
            // Validity: QueryRecord contains uppercase bases
            && let Some(r) = query.nucleotides.find_substring(b"ATG")
            && query.nucleotides.len() - r.start >= ref_id_data.length
        {
            (r.start, query.nucleotides.slice(r.start..))
        } else {
            (0, query.nucleotides.as_view())
        }
    }

    /// If the `repairable_end_limit` rule is enabled with a non-zero limit $L$,
    /// then extend the ends of the alignment if they are within the limit.
    ///
    /// Specifically, if at most $L$ bases were clipped from both sequences on
    /// the left, then extend the alignment to include these bases. Similarly,
    /// if at most $L$ bases were clipped from both sequences on the right, then
    /// extend the alignment as well.
    ///
    /// The score of the alignment is not altered.
    fn rule_repairable_ends(&self, genome_aln: &mut Alignment<u32>) {
        if let Some(limit) = self.data.rules.repairable_end_limit {
            let unaligned_pre = genome_aln.query_range.start.min(genome_aln.ref_range.start);
            if genome_aln.ref_range.start <= limit {
                genome_aln.extend_left(unaligned_pre);
            }

            let unaligned_ref_post = genome_aln.uanligned_ref_tail();
            let unaligned_post = genome_aln.unaligned_query_tail().min(unaligned_ref_post);
            if unaligned_ref_post <= limit {
                genome_aln.extend_right(unaligned_post);
            }
        }
    }
}

/// An extension trait for [`Alignment`] offering methods specific to
/// DAIS-ribosome.
trait AlignmentExt {
    /// Extends the alignment to the left by `by_length` bases.
    ///
    /// In the case where clipping at the start of the query and reference
    /// sequences both occurred, the left side of the alignment can be
    /// "extended" by converting some number of these bases to match states.
    /// This method does not alter the score field.
    fn extend_left(&mut self, by_length: usize);

    /// Extends the alignment to the right by `by_length` bases.
    ///
    /// In the case where clipping at the end of the query and reference
    /// sequences both occurred, the right side of the alignment can be
    /// "extended" by converting some number of these bases to match states.
    /// This method does not alter the score field.
    fn extend_right(&mut self, by_length: usize);

    /// Returns the number of unaligned bases at the end of the query.
    fn unaligned_query_tail(&self) -> usize;

    /// Returns the number of unaligned bases at the end of the reference.
    fn uanligned_ref_tail(&self) -> usize;

    /// Returns an owned copy of the aligned query range.
    fn aln_query_range(&self) -> Range<usize>;
}

impl<T> AlignmentExt for Alignment<T> {
    fn aln_query_range(&self) -> Range<usize> {
        self.query_range.clone()
    }

    fn unaligned_query_tail(&self) -> usize {
        self.query_len - self.query_range.end
    }

    fn uanligned_ref_tail(&self) -> usize {
        self.ref_len - self.ref_range.end
    }

    fn extend_left(&mut self, by_length: usize) {
        self.query_range.start -= by_length;
        self.ref_range.start -= by_length;
        self.states.extend_left(by_length);
    }

    fn extend_right(&mut self, by_length: usize) {
        self.query_range.end += by_length;
        self.ref_range.end += by_length;
        self.states.extend_right(by_length);
    }
}

impl<'a> Product<'a> {
    /// Returns true if a required start codon is required for the exons, a
    /// long-enough match state occurs in the [`Product`] to span this codon,
    /// and yet the codon is not equal to the required one.
    ///
    /// If any of these conditions are not true, then `false` is returned.
    ///
    /// `T` and `U` are treated equivalently for the purposes of this check.
    ///
    /// ## Validity
    ///
    /// The `query` should be the same query which the alignment used to create
    /// `self` was formed from.
    pub(crate) fn missing_required_start(&self, query: &QueryRecord) -> bool {
        if let Some(required) = self.product_spec.exons.required_start
            // Note that the first product range is either a match or deletion
            && let Some(CdsStateRange::M(m)) = self.product_ranges.first()
            && m.cds_range.len() >= 3
        {
            // Validity: query_range always refers to valid indices in query. At
            // least 3 residues exist since cds_range and query_range are the
            // same length
            let mut first: [u8; 3] = *query.nucleotides[m.query_range.start..].first_chunk().expect("The length of the query_range should be at least 3, and the query_range should not refer to out of bounds indices");

            // Convert U to T for the purpose of identifying the required start
            first = first.map(|b| b.recode_base(RecodeDNAStrat::AnyToAcgtnNoGapsUpper));

            first != required
        } else {
            false
        }
    }

    /// Condenses adjacent deletions in the coding sequence.
    ///
    /// If a deletion spans multiple exons, the intersection algorithm will
    /// split it into separate ranges. This method combines them again for
    /// downstream correctness.
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

    /// Fixes indel frames by repositioning insertions and deletions to in-frame
    /// boundaries.
    ///
    /// This corrects in-frame indels (length divisible by 3) that occur at
    /// out-of-frame positions by shifting them left or right based on codon
    /// usage statistics. The algorithm mirrors the logic in
    /// `codonCorrectStats.pl`, albeit for coordinate math instead of strings.
    pub(crate) fn fix_frames(&mut self, query: &QueryRecord) {
        let len = self.product_ranges.len();

        // We iterate by index because we need to access neighbors
        let mut i = 0;
        while i < len {
            match &self.product_ranges[i] {
                CdsStateRange::I(ins) => {
                    let codon_shift = ins.cds_index.codon_shift();

                    // Only correct in-frame insertions at out-of-frame positions
                    if codon_shift != 0 && ins.len() % 3 == 0 {
                        self.fix_insertion_frame(i, query);
                    }
                }
                CdsStateRange::D(del) => {
                    // The inclusive 0-based index is the upstream 1-based position
                    let frame = del.cds_range.start % 3;

                    // Only correct in-frame deletions at out-of-frame positions
                    if frame != 0 && del.len() % 3 == 0 {
                        self.fix_deletion_frame(i, query);
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
    fn fix_insertion_frame(&mut self, idx: usize, query: &QueryRecord) -> Option<()> {
        let Product {
            product_ranges,
            product_spec,
            ..
        } = self;
        let query = &query.nucleotides;

        // Insertions cannot be split over exons, and in practice many scoring
        // schemes tend to avoid adjacent insertions and deletions, so in
        // practice left_match and right_match occur at idx-1 and idx+1
        // respectively. But in case either is a deletion, we perform a search.
        let Some((left_match, CdsStateRange::I(ins), right_match)) = Self::partition_states(product_ranges, idx) else {
            return None;
        };

        let codon_shift = ins.cds_index.codon_shift();
        let codon_position = ins.cds_index.to_aa_idx().right_pos();
        let insert_len = ins.len();

        // for weight lookup
        let insert_seq = &query[ins.query_range.clone()];

        if codon_shift == 2 {
            // A2: insertion after 2nd codon base
            //
            // A2/L2: shift insert left 2
            //   New codon: last 2 bases of insert + cp3
            //
            // A2/R1: shift insert right 1
            //   New codon: cp1 + cp2 + first 1 base of insert
            let cp1 = *query.get(left_match.query_range.end - 2)?;
            let cp2 = *query.get(left_match.query_range.end - 1)?;
            let cp3 = *query.get(right_match.query_range.start)?;

            let a2l2 = build_discontiguous_codon(insert_seq[insert_len - 2], insert_seq[insert_len - 1], cp3);
            let a2r1 = build_discontiguous_codon(cp1, cp2, insert_seq[0]);

            // Validity: The codons are uppercase because they are derived from
            // bases in query
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
        } else if codon_shift == 1 {
            // A1 insertion: insertion after 1st base of codon
            //
            // A1/L1: shift insert left 1
            //   New codon: last 1 base of insert + cp2 + cp3
            //
            // A1/R2: shift insert right 2
            //   New codon: cp1 + first 2 bases of insert

            let cp1 = *query.get(left_match.query_range.end - 1)?;
            let cp2 = *query.get(right_match.query_range.start)?;
            let cp3 = *query.get(right_match.query_range.start + 1)?;

            let a1l1 = build_discontiguous_codon(insert_seq[insert_len - 1], cp2, cp3);
            let a1r2 = build_discontiguous_codon(cp1, insert_seq[0], insert_seq[1]);

            // Validity: build_discontiguous_codon ensures validity requirements
            // are met
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
    fn fix_deletion_frame(&mut self, idx: usize, query: &QueryRecord) -> Option<()> {
        let Product {
            product_ranges,
            product_spec,
            ..
        } = self;
        let query = &query.nucleotides;

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
            let pivot = build_discontiguous_codon(
                *query.get(left_match.query_range.end - 1)?,
                *query.get(right_match.query_range.start)?,
                *query.get(right_match.query_range.start + 1)?,
            );

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
            let pivot = build_discontiguous_codon(
                *query.get(left_match.query_range.end - 2)?,
                *query.get(left_match.query_range.end - 1)?,
                *query.get(right_match.query_range.start)?,
            );

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
}

/// Combines three discontinuous bases into a codon, automatically converting
/// `U` to `T`.
///
/// No other bytes are altered, and case is not changed.
fn build_discontiguous_codon(b1: u8, b2: u8, b3: u8) -> [u8; 3] {
    [b1, b2, b3].map(|base| if base == b'U' { b'T' } else { base })
}

/// An extension trait for slices, providing functionality specific to
/// DAIS-ribosome.
trait SliceExt<T> {
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
