use crate::{
    config::{ProductSpec, toml::Formatting},
    data::{
        QueryRecord,
        ranges::{CdsStateRange, InsertionIdx, StateRange},
    },
    hashing::nt_id_iupac,
    outputs::{ComputedGenome, ComputedGenomeInsertion},
};
use std::ops::Range;
use zoe::prelude::*;

/// The genome alignments and products for a single query against all reference
/// IDs.
///
/// The genome information can be materialized into [`ComputedGenome`] with
/// [`GenomeAndProductStates::materialize_genome`]. Each product can be
/// materialized into [`ComputedProduct`] with [`Product::materialize`].
///
/// [`ComputedProduct`]: crate::outputs::ComputedProduct
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

/// The genome alignment and products for a single query against a single
/// reference.
///
/// The genome information can be materialized into [`ComputedGenome`] with
/// [`materialize_genome`]. Each product can be materialized into
/// [`ComputedProduct`] with [`Product::materialize`].
///
/// Many [`GenomeAndProductStates`] are stored in [`RibosomeOutput`], which
/// holds the alignments for all reference IDs.
///
/// [`materialize_genome`]: GenomeAndProductStates::materialize_genome
/// [`ComputedProduct`]: crate::outputs::ComputedProduct
#[derive(Debug)]
pub struct GenomeAndProductStates<'a> {
    /// The ID for the reference group which was aligned against.
    pub reference_id: &'a str,

    /// The length of the reference sequence.
    pub ref_len: usize,

    /// Genome alignment to nucleotide reference sequence expressed as [`StateRange`]
    pub genome_aln_states: Vec<StateRange>,

    /// The range of the stop extension within the query, if present.
    pub stop_extension_query_range: Option<Range<usize>>,

    /// The number of bases in the reference sequence that were not aligned
    /// against in the beginning (i.e., not included in `genome_aln_states`).
    pub leading_ref_unaligned: usize,

    /// The number of bases in the reference sequence that were not aligned
    /// against at the end (i.e., not included in `genome_aln_states`).
    pub trailing_ref_unaligned: usize,

    /// Contains all relevant product data, including the protein name.
    pub products: Vec<Product<'a>>,
}

impl<'a> GenomeAndProductStates<'a> {
    /// Computes the output data for this genome, materializing all ranges into
    /// sequences using `query`.
    pub fn materialize_genome(&self, query: &QueryRecord) -> ComputedGenome {
        let query = query.nucleotides();

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

                    // Validity: slice is from QueryRecord
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

            // Validity: slice is from QueryRecord
            insertions.push(ComputedGenomeInsertion::new(nt_insertion_idx, slice));
        }

        // Validity: genome_seq contains unaligned uppercase IUPAC since
        // QueryRecord does
        let genome_id = nt_id_iupac(&genome_seq);

        let genome_length = genome_seq.len();

        ComputedGenome {
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
/// This can be materialized into [`ComputedProduct`] using [`materialize`].
///
/// Many [`Product`] values are stored in [`GenomeAndProductStates`], which
/// holds the alignments for all the protein products. Many of these are stored
/// in [`RibosomeOutput`], which holds the alignments for all reference IDs.
///
/// [`materialize`]: Product::materialize
/// [`ComputedProduct`]: crate::outputs::ComputedProduct
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
    /// exons. The contained ranges will have non-zero length.
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
    ///
    /// If `Some`, the range will have length at least 3. The last three indices
    /// will correspond to the stop codon in `query`. This can only be set if
    /// `product_ranges` ends in a match state.
    pub stop_extension_query_range: Option<Range<usize>>,
}
