use crate::{
    config::Formatting,
    data::{
        ComputedGenomeInsertion, PrecomputedGenomeData, QueryRecord,
        products::Product,
        ranges::{InsertionIdx, StateRange},
    },
    hashing::nt_id,
};
use std::ops::Range;
use zoe::prelude::*;

#[derive(Debug)]
pub struct RibosomeOutput<'a> {
    /// Original query record
    pub query:             QueryRecord,
    /// Both genome and protein product alignment states
    pub(crate) states:     Vec<GenomeAndProductStates<'a>>,
    /// Output formatting rules, parsed from the TOML configuration.
    pub(crate) formatting: &'a Formatting,
}

#[derive(Debug)]
pub(crate) struct GenomeAndProductStates<'a> {
    /// Reference ID
    pub(crate) reference_id: &'a str,

    /// The length of the reference sequence.
    pub(crate) ref_len: usize,

    /// Genome alignment to nucleotide reference sequence expressed as [`StateRange`]
    pub(crate) genome_aln_states: Vec<StateRange>,

    /// The range of the stop extension within the query, if present.
    pub(crate) stop_extension_query_range: Option<Range<usize>>,

    /// The number of bases in the reference sequence that were not aligned
    /// against in the beginning (i.e., not included in `genome_aln_states`).
    pub(crate) leading_ref_unaligned: usize,

    /// The number of bases in the reference sequence that were not aligned
    /// against at the end (i.e., not included in `genome_aln_states`).
    pub(crate) trailing_ref_unaligned: usize,

    /// Contains all relevant product data, including the protein name.
    pub(crate) products: Vec<Product<'a>>,
}

impl<'a> GenomeAndProductStates<'a> {
    /// Lazily compute and cache genome data from genome alignment states.
    pub fn materialize_genome(&self, query: &Nucleotides) -> PrecomputedGenomeData {
        let mut genome_seq = Nucleotides::new();
        let mut genome_aln = Nucleotides::from(vec![b'.'; self.leading_ref_unaligned]);
        let mut insertions = Vec::new();
        let mut has_insertion = false;

        for state in &self.genome_aln_states {
            match state {
                StateRange::M(m) => {
                    let slice = &query[m.query_range.clone()];
                    genome_seq.extend_from_slice(slice);
                    genome_aln.extend_from_slice(slice);
                }
                StateRange::I(ins) => {
                    let slice = &query[ins.query_range.clone()];
                    genome_seq.extend_from_slice(slice);
                    has_insertion = true;
                    insertions.push(ComputedGenomeInsertion::new(ins.ref_index, slice));
                }
                StateRange::D(del) => {
                    genome_aln.pad_end(b'-', del.ref_range.len());
                }
            }
        }

        if let Some(ref ext_range) = self.stop_extension_query_range {
            let nt_insertion_idx = InsertionIdx::from_right_idx(self.ref_len);
            let slice = &query[ext_range.clone()];
            genome_seq.extend_from_slice(slice);
            has_insertion = true;
            insertions.push(ComputedGenomeInsertion::new(nt_insertion_idx, slice));
        }

        let genome_id = nt_id(&genome_seq);
        let genome_length = genome_seq.len();

        PrecomputedGenomeData {
            genome_id,
            genome_length,
            has_insertion,
            genome_seq,
            genome_aln,
            insertions,
            trailing_ref_unaligned: self.trailing_ref_unaligned,
        }
    }
}
