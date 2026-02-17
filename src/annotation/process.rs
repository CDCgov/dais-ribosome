use std::ops::Range;
use std::sync::OnceLock;
use zoe::{alignment::Alignment, data::types::nucleotides::CodonExtension, prelude::*};

use crate::{
    annotation::error::RibosomeError,
    data::{
        GenomeAndProductStates, RibosomeOutput,
        ctype::ReferenceGroup,
        query::QueryRecord,
        ranges::{InsertionRange, StateRange},
    },
};

use super::*;

impl<'a> AnnotationModule<'a> {
    pub fn process(&self, query: QueryRecord) -> Result<RibosomeOutput<'_>, RibosomeError> {
        let Some(reference_data) = self.ctype_map.get(&query.ctype) else {
            return Err(RibosomeError::UnimplementedCtype(query.ctype.to_string()));
        };

        let mut states = Vec::with_capacity(reference_data.len());

        for ref_id_data in reference_data.iter() {
            let (query_ori_offset, query_seq) = self.rule_chew_to_start(&query, ref_id_data);

            let Some(mut genome_aln) = ref_id_data.best_alignment(&query_seq) else {
                return Err(RibosomeError::Unmappable(query.id.to_string()));
            };

            //eprintln!("{}\t{}", query.id, genome_aln.states);

            self.rule_repairable_ends(&mut genome_aln);

            let stop_extension = self.rule_stop_extension(&query_seq, &genome_aln);
            let mut genome_aln_states = StateRange::state_ranges_from_aligment(&genome_aln);
            let mut products = Vec::with_capacity(ref_id_data.proteins.len());

            for product in ref_id_data.iter_proteins() {
                let mut product_ranges = product.make_product_ranges(&genome_aln_states);
                if product_ranges.missing_required_start(query_seq) {
                    continue;
                }

                product_ranges.condense_deletions();
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
                    match state {
                        StateRange::M(m) => m.shift_query_right(query_ori_offset),
                        StateRange::I(i) => i.shift_query_right(query_ori_offset),
                        _ => {}
                    }
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

    fn rule_stop_extension<'b>(
        &self, query_seq: &'b NucleotidesView<'b>, genome_aln: &Alignment<u32>,
    ) -> Option<InsertionRange> {
        if self.data.rules.list_contig_stop_extension
            && genome_aln.unaligned_query_tail() >= 3
            && let Some(last_codon) = query_seq.slice(genome_aln.aln_query_range()).get_tail_codon()
            && last_codon.is_amino_acid()
            && let Some(stop_codon_index) = query_seq.slice(genome_aln.query_range.end..).find_next_aa_in_frame(b'*')
        {
            // For half-open 0-based ranges, the exclusive end is the inclusive starting index
            let start_index = genome_aln.query_range.end;
            Some(InsertionRange {
                upstream_ref_index: genome_aln.ref_range.end - 1,
                query_range:        start_index..start_index + stop_codon_index + 3,
            })
        } else {
            None
        }
    }

    fn rule_chew_to_start<'b>(
        &self, query: &'b QueryRecord, ref_id_data: &ReferenceGroup<'_>,
    ) -> (usize, NucleotidesView<'b>) {
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

trait AlignmentExt {
    fn extend_left(&mut self, by_length: usize);
    fn extend_right(&mut self, by_length: usize);
    fn unaligned_query_tail(&self) -> usize;
    fn uanligned_ref_tail(&self) -> usize;
    fn aln_query_range(&self) -> Range<usize>;
}

impl<T> AlignmentExt for zoe::alignment::Alignment<T> {
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
        if by_length > 0 {
            self.query_range.start -= by_length;
            self.ref_range.start -= by_length;
            self.states.prepend_inc_op(by_length, b'M');
        }
    }

    fn extend_right(&mut self, by_length: usize) {
        if by_length > 0 {
            self.query_range.end += by_length;
            self.ref_range.end += by_length;
            self.states.add_inc_op(by_length, b'M');
        }
    }
}
