use crate::{
    annotation::{AnnotationModule, error::RibosomeError},
    data::{
        GenomeAndProductStates, RibosomeOutput,
        ctype::ReferenceGroup,
        query::QueryRecord,
        ranges::{InsertionRange, StateRange},
    },
};
use std::{ops::Range, sync::OnceLock};
use zoe::{alignment::Alignment, data::types::nucleotides::CodonExtension, prelude::*};

impl<'a> AnnotationModule<'a> {
    /// Processes a single query TODO: what does it return?
    pub fn process(&self, query: QueryRecord) -> Result<RibosomeOutput<'_>, RibosomeError> {
        // Get the corresponding reference information for the compound type of
        // the query
        let Some(reference_data) = self.ctype_map.get(&query.ctype) else {
            return Err(RibosomeError::UnimplementedCtype(query.ctype.to_string()));
        };

        let mut states = Vec::with_capacity(reference_data.len());

        for ref_id_data in reference_data.iter() {
            let (query_ori_offset, query_seq) = self.rule_chew_to_start(&query, ref_id_data);

            // Get the alignment to the best reference
            let Some(mut genome_aln) = ref_id_data.best_alignment(&query_seq) else {
                return Err(RibosomeError::Unmappable(query.id.to_string()));
            };

            //eprintln!("{}\t{}", query.id, genome_aln.states);

            // Extend the left and right side of the alignments
            self.rule_repairable_ends(&mut genome_aln);

            let stop_extension = self.rule_stop_extension(&query_seq, &genome_aln);
            let mut genome_aln_states = StateRange::state_ranges_from_aligment(&genome_aln);
            let mut products = Vec::with_capacity(ref_id_data.proteins.len());

            for product in &ref_id_data.proteins {
                let mut product_ranges = product.make_product_ranges(&genome_aln_states);
                if product_ranges.missing_required_start(query_seq) {
                    continue;
                }

                product_ranges.condense_deletions();
                // Validity: QueryRecord contains unaligned, uppercase IUPAC
                // bases
                product_ranges.fix_frames(query_seq);
                product_ranges.add_query_coords(query_ori_offset);

                products.push(product_ranges);
            }

            // Push stop extension into every product whose last exon ends at
            // the extension's reference position (matching Perl $pMax == $max).
            if let Some(mut ext) = stop_extension {
                if query_ori_offset > 0 {
                    ext.shift_query_right(query_ori_offset);
                }

                let ext_ref_end = ext.upstream_ref_index + 1;
                for product in &mut products {
                    if let Some(last_exon) = product.product_spec.exons.coords.last()
                        && last_exon.ref_range.end == ext_ref_end
                    {
                        product.stop_extension_query_range = Some(ext.query_range.clone());
                    }
                }
            }

            if query_ori_offset > 0 {
                for state in &mut genome_aln_states {
                    state.shift_query_right(query_ori_offset);
                }
            }

            states.push(GenomeAndProductStates {
                reference_id: &ref_id_data.reference_id,
                ref_len: ref_id_data.length,
                genome_aln_states,
                products,
                computed_genome: OnceLock::new(),
            });
        }

        Ok(RibosomeOutput {
            query,
            states,
            formatting: &self.data.formatting,
        })
    }

    // TODO: Slicing the query and then searching in frame could cause issues if
    // the slice disrupts the frame

    fn rule_stop_extension<'b>(
        &self, query_seq: &'b NucleotidesView<'b>, genome_aln: &Alignment<u32>,
    ) -> Option<InsertionRange> {
        if self.data.rules.list_contig_stop_extension
            && genome_aln.unaligned_query_tail() >= 3
            && let Some(last_aligned_codon) = query_seq.slice(genome_aln.aln_query_range()).get_tail_codon()
            && last_aligned_codon.is_amino_acid()
            && let Some(stop_codon_index) = query_seq.slice(genome_aln.query_range.end..).find_next_aa_in_frame(b'*')
        {
            // The exclusive end of the last alignment range is the inclusive
            // start of the insertion
            let start_index = genome_aln.query_range.end;

            // TODO: WE HAVE ISSUES!
            let end_index = start_index + stop_codon_index + 3;

            Some(InsertionRange {
                // The last aligned residue is ref_range.end - 1 (converting
                // exclusive to inclusive). The insertion occurs after this
                // residue.
                upstream_ref_index: genome_aln.ref_range.end - 1,
                query_range:        start_index..start_index + stop_codon_index + 3,
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
        // Validity: QueryRecord guarantees U has been replaced with T, so ATG
        // is the only possible start codon

        if self.data.rules.chew_to_start
            && query.nucleotides.len() > ref_id_data.length
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
        self.states.prepend_inc_op(by_length, b'M');
    }

    fn extend_right(&mut self, by_length: usize) {
        self.query_range.end += by_length;
        self.ref_range.end += by_length;
        self.states.add_inc_op(by_length, b'M');
    }
}
