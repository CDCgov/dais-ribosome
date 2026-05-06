use crate::{
    config::{product_spec::ProductSpec, toml::Formatting},
    data::{
        QueryRecord,
        ranges::{CdsStateRange, InsertionIdx, StateRange},
    },
    hashing::nt_id,
    outputs::{ComputedGenomeInsertion, PrecomputedGenomeData},
};
use std::ops::Range;
use zoe::prelude::*;

#[derive(Debug)]
pub struct RibosomeOutput<'a> {
    /// Original query record
    pub query:          QueryRecord,
    /// Both genome and protein product alignment states
    pub states:         Vec<GenomeAndProductStates<'a>>,
    /// A vector of any reference IDs which failed to be aligned against
    pub failed_ref_ids: Vec<String>,
    /// Output formatting rules, parsed from the TOML configuration.
    pub formatting:     &'a Formatting,
}

#[derive(Debug)]
pub struct GenomeAndProductStates<'a> {
    /// Reference ID
    pub reference_id: &'a str,

    /// The length of the reference sequence.
    pub(crate) ref_len: usize,

    /// Genome alignment to nucleotide reference sequence expressed as [`StateRange`]
    pub genome_aln_states: Vec<StateRange>,

    /// The range of the stop extension within the query, if present.
    pub(crate) stop_extension_query_range: Option<Range<usize>>,

    /// The number of bases in the reference sequence that were not aligned
    /// against in the beginning (i.e., not included in `genome_aln_states`).
    pub(crate) leading_ref_unaligned: usize,

    /// The number of bases in the reference sequence that were not aligned
    /// against at the end (i.e., not included in `genome_aln_states`).
    pub(crate) trailing_ref_unaligned: usize,

    /// Contains all relevant product data, including the protein name.
    pub products: Vec<Product<'a>>,
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

/// The aligned ranges for a single query against a single reference, using the
/// exons for one of the protein products.
///
/// Many [`Product`] values are stored in [`GenomeAndProductStates`], which
/// holds the alignments for all the protein products. Many of these are stored
/// in [`RibosomeOutput`], which holds the alignments for all reference IDs.
/// Many [`RibosomeOutput`] are generated and written in `main.rs` for each
/// query.
///
/// [`GenomeAndProductStates`]: crate::outputs::GenomeAndProductStates
/// [`RibosomeOutput`]: crate::outputs::RibosomeOutput
#[derive(Debug)]
pub struct Product<'a> {
    /// The information for the protein product being aligned against, including
    /// the name and exons.
    pub(crate) product_spec: &'a ProductSpec,

    /// The ranges within the exons that the query covers. This is initially
    /// formed by intersecting the query ranges with the exon ranges, then is
    /// edited.
    ///
    /// The [`CdsStateRange`] values are ordered and partition the
    /// aligned-against range of the coding sequence. Overlapping exons can
    /// cause the query ranges to have repeated indices and not be in order, so
    /// no guarantees can be made for those fields.
    ///
    /// This field does not begin or end with [`CdsStateRange::I`]. This field
    /// may be empty if there is no intersection between the query and the
    /// exons.
    pub product_ranges: Vec<CdsStateRange>,

    /// The number of bases in the coding sequence that were not aligned against
    /// in the beginning (i.e., not included in `product_ranges`).
    ///
    /// If `product_ranges` is empty, then this field is 0, and the unaligned
    /// bases are counted in `trailing_cds_unaligned`.
    pub leading_cds_unaligned: usize,

    /// The number of bases in the coding sequence that were not aligned against
    /// at the end (i.e., not included in `product_ranges`).
    ///
    /// Materialization may use an _increased_ version of this if it truncates
    /// the alignment at the first stop codon encountered.
    pub trailing_cds_unaligned: usize,

    /// If this product's last exon ends at the stop extension position, this
    /// holds the query range of the stop extension nucleotides.
    pub stop_extension_query_range: Option<Range<usize>>,
}
