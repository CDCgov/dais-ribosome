//! The core process of DAIS-ribosome for creating a range-based
//! [`RibosomeOutput`] annotation result from an input [`QueryRecord`].

use crate::{
    AlignmentStatesExt,
    annotation::intersection,
    config::annotation_module::{AnnotationModule, ReferenceGroup},
    data::{
        QueryRecord,
        ranges::{CdsStateRange, InsertionIdx, InsertionRange, RangeExt, StateRange},
    },
    errors::RibosomeError,
    outputs::{GenomeAndProductStates, Product, RibosomeOutput},
    ranges::CdsInsertionRange,
};
use std::ops::Range;
use zoe::{alignment::Alignment, data::types::nucleotides::CodonExtension, prelude::*};

impl<'a> AnnotationModule<'a> {
    /// Processes a single query, returning [`RibosomeOutput`] containing all
    /// the genome alignments against the relevant references, as well as the
    /// protein products formed for each reference.
    pub fn process(&self, query: QueryRecord) -> Result<RibosomeOutput<'_>, RibosomeError> {
        // Get the corresponding reference information for the compound type of
        // the query
        let Some(reference_data) = self.ctype_map.get(query.ctype()) else {
            return Err(RibosomeError::UnimplementedCtype(query.into_ctype().into()));
        };

        let mut states = Vec::with_capacity(reference_data.len());

        let mut failed_ref_ids = Vec::new();

        for ref_id_data in reference_data.iter() {
            let (query_ori_offset, chewed_query) = self.rule_chew_to_start(&query, ref_id_data);

            // Get the alignment to the best reference
            let Some(mut genome_aln) = self.best_alignment(ref_id_data, &chewed_query) else {
                failed_ref_ids.push(ref_id_data.reference_id.clone());
                continue;
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

            let mut products = Vec::with_capacity(ref_id_data.product_specs.len());

            for product_spec in &ref_id_data.product_specs {
                // Validity: requirements met based on
                // state_ranges_from_aligment guarantees
                let mut product = intersection::form_product(&genome_aln_states, product_spec);

                // Shift indels to fix their frames. Validity: this is called
                // before condense_deletions.
                product.fix_frames(&query);

                // Condense any remaining deletions
                product.condense_deletions();

                // Validity: the same `query_seq` is passed as was used to form
                // `genome_aln_states`
                if product.missing_required_start(&query) {
                    continue;
                }

                // Add the stop extension if applicable
                if let Some(ext) = &stop_extension
                    && product.product_spec.exons.last().ref_range.end == ext.ref_index.right()
                {
                    // stop_extension is only Some if the genome alignment
                    // extends to the end of the reference. Above we ensure the
                    // exons extend to the end of the reference. Hence, the
                    // product alignment (intersection of the two) will extend
                    // to the end of the reference.
                    debug_assert_eq!(product.rpad, 0);

                    let stop_extension = CdsInsertionRange {
                        cds_index:   InsertionIdx::from_right_idx(product.product_spec.exons.cds_len()),
                        query_range: ext.query_range.clone(),
                    };

                    product.product_ranges.push(CdsStateRange::I(stop_extension));
                }

                products.push(product);
            }

            let ref_len = ref_id_data.length;

            let lpad = genome_aln_states
                .iter()
                .find_map(|s| match s {
                    StateRange::M(m) => Some(m.ref_range.start),
                    StateRange::D(d) => Some(d.ref_range.start),
                    _ => None,
                })
                .unwrap_or(0);

            let rpad = ref_len
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
                lpad,
                rpad,
                products,
                stop_extension_query_range: stop_extension.map(|ins| ins.query_range),
            });
        }

        Ok(RibosomeOutput {
            query,
            states,
            failed_ref_ids,
            formatting: self.formatting,
        })
    }

    /// Computes the stop extension if the [`list_contig_stop_extension`] rule
    /// is set.
    ///
    /// If `Some`, the returned range will have length at least 3. The last
    /// three indices will correspond to the stop codon in `query`.
    ///
    /// [`list_contig_stop_extension`]:
    ///     crate::config::toml::Rules::list_contig_stop_extension
    fn rule_stop_extension(&self, query: &QueryRecord, genome_aln: &Alignment<u32>) -> Option<InsertionRange> {
        let query_seq = query.nucleotides();

        if self.rules.list_contig_stop_extension
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

    /// Applies the [`chew_to_start`] rule if it is set.
    ///
    /// The first value returned is the starting index of the returned slice
    /// within the query. If the rule cannot be applied or is not set, then the
    /// full query is returned and with a starting index of 0.
    ///
    /// [`chew_to_start`]: crate::config::toml::Rules::chew_to_start
    fn rule_chew_to_start<'b>(
        &self, query: &'b QueryRecord, ref_id_data: &ReferenceGroup<'_>,
    ) -> (usize, NucleotidesView<'b>) {
        let query = query.nucleotides();

        if self.rules.chew_to_start
            && query.len() > ref_id_data.length
            // Validity: QueryRecord contains uppercase bases
            && let Some(r) = query.find_substring(b"ATG")
            && query.len() - r.start >= ref_id_data.length
        {
            (r.start, query.slice(r.start..))
        } else {
            (0, query.as_view())
        }
    }

    /// Applies the repairable ends rule with the given limit in
    /// [`repairable_end_limit`].
    ///
    /// The score of the alignment is not altered. The alignment may have match
    /// states added and soft clipping removed, but no other changes will occur.
    ///
    /// [`repairable_end_limit`]:
    ///     crate::config::toml::Rules::repairable_end_limit
    fn rule_repairable_ends(&self, genome_aln: &mut Alignment<u32>) {
        let limit = self.rules.repairable_end_limit;

        if limit > 0 {
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
        let Some(required) = self.product_spec.exons.required_start else {
            // The specs do not require a start codon, so return false
            return false;
        };

        // Get the index of the start codon in query coordinates, if available
        let start_codon_idx_in_query = match self.product_ranges.first() {
            Some(CdsStateRange::M(state)) => (state.cds_range.start == 0).then_some(state.query_range.start),
            // Leading insertions shouldn't be possible, but we can handle them anyways
            Some(CdsStateRange::I(state)) => state.cds_index.at_start().then_some(state.query_range.start),
            // Deletion or empty ranges implies start codon was partly deleted
            // or clipped
            Some(CdsStateRange::D(_)) | None => None,
        };

        let Some(start_codon_idx_in_query) = start_codon_idx_in_query else {
            return true;
        };

        // We use get_slice since technically an empty match state could cause
        // query_range.start to be out of bounds, but this will never happen
        let Some(mut first_cds_codon) = query
            .nucleotides()
            .get_slice(start_codon_idx_in_query..)
            .and_then(|slice| slice.get_first_codon())
        else {
            // The entire CDS is less than a single codon in length, so the required start codon is missing
            return true;
        };

        // Convert U to T for the purpose of identifying the required start
        first_cds_codon = first_cds_codon.map(|b| if b == b'U' { b'T' } else { b });

        // Check for equality
        first_cds_codon != required
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
                // Validity: product_ranges partitions the aligned-against
                // coding sequence, so adjacent deletions in product_ranges are
                // also adjacent in the coding sequence
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
}
