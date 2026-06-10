use crate::{
    config::{annotation_module::ReferenceGroup, toml::Formatting},
    data::{
        QueryRecord,
        ranges::{CdsStateRange, InsertionIdx, StateRange},
    },
    hashing::nt_id_iupac,
    outputs::{ComputedGenome, ComputedGenomeInsertion},
    ranges::InsertionRange,
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

    /// The genome alignment to the nucleotide reference sequence.
    ///
    /// This will begin and end with match states.
    pub genome_aln_states: Vec<StateRange>,

    /// The range of the stop extension within the query, if present.
    pub stop_extension_query_range: Option<Range<usize>>,

    /// The amount of left padding in the genome alignment.
    ///
    /// This is the number of bases in the reference sequence that were not
    /// aligned against at the beginning (i.e., not included in
    /// `genome_aln_states`).
    pub lpad: usize,

    /// The amount of right padding in the genome alignment.
    ///
    /// This is the number of bases in the reference sequence that were not
    /// aligned against at the end (i.e., not included in `genome_aln_states`).
    pub rpad: usize,

    /// Contains all relevant product data, including the protein name.
    pub products: Vec<Product<'a>>,
}

impl<'a> GenomeAndProductStates<'a> {
    /// Creates a new [`GenomeAndProductStates`] from the given reference
    /// information, states, stop extension, and products. The stop extension
    /// must already be added to the products if applicable.
    ///
    /// ## Panics
    ///
    /// In debug mode, this panics if `genome_aln_states` does not begin and end
    /// with a match state.
    pub(crate) fn new(
        references: &'a ReferenceGroup, genome_aln_states: Vec<StateRange>, stop_extension: Option<InsertionRange>,
        products: Vec<Product<'a>>,
    ) -> Self {
        #[cfg(debug_assertions)]
        {
            let (Some(first), Some(last)) = (genome_aln_states.first(), genome_aln_states.last()) else {
                panic!("genome_aln_states cannot be empty");
            };

            if !first.is_match() {
                panic!("genome_aln_states must start with a match state");
            }

            if !last.is_match() {
                panic!("genome_aln_states must end with a match state");
            }
        }

        let ref_len = references.length;

        let lpad = genome_aln_states
            .first()
            .map_or(0, |first_state| first_state.begin_ref_coord());
        let rpad = ref_len - genome_aln_states.last().map_or(0, |last_state| last_state.end_ref_coord());

        Self {
            reference_id: &references.reference_id,
            ref_len,
            genome_aln_states,
            stop_extension_query_range: stop_extension.map(|ins| ins.query_range),
            lpad,
            rpad,
            products,
        }
    }

    /// Computes the output data for this genome, materializing all ranges into
    /// sequences using `query`.
    pub fn materialize_genome(&self, query: &QueryRecord) -> ComputedGenome {
        let query = query.nucleotides();

        let mut genome_seq = Nucleotides::new();
        let mut genome_aln = Nucleotides::from(vec![b'.'; self.lpad]);
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
        // QueryRecord does. It is non-empty since genome_aln_states begins with
        // a match state.
        let genome_id = nt_id_iupac(&genome_seq);

        let genome_length = genome_seq.len();

        ComputedGenome {
            genome_id,
            genome_length,
            has_insertion,
            genome_seq,
            genome_aln,
            insertions,
            genome_aln_rpad: self.rpad,
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
    pub name: &'a str,

    /// The ranges within the exons that the query covers. This is initially
    /// formed by intersecting the query ranges with the exon ranges, then is
    /// edited.
    ///
    /// The [`CdsStateRange`] values are ordered and partition the
    /// aligned-against range of the coding sequence. Overlapping exons can
    /// cause the query ranges to have repeated indices and not be in order, so
    /// no guarantees can be made for those fields.
    ///
    /// This may end in a trailing insertion, which represents a stop extension.
    ///
    /// This field may be empty if there is no intersection between the query
    /// and the exons. The contained ranges will have non-zero length.
    pub product_ranges: Vec<CdsStateRange>,

    /// The amount of left padding in CDS coordinates for the product.
    ///
    /// This is the number of bases in the coding sequence that were not aligned
    /// against in the beginning (i.e., not included in `product_ranges`).
    ///
    /// If `product_ranges` is empty, then this field is 0, and the unaligned
    /// bases are counted in `rpad`.
    pub lpad: usize,

    /// The number of bases in the coding sequence that were aligned against
    /// (i.e., included in `product_ranges`).
    pub cds_aligned_len: usize,

    /// The amount of right padding in CDS coordinates for the product.
    ///
    /// This is the number of bases in the coding sequence that were not aligned
    /// against at the end (i.e., not included in `product_ranges`).
    ///
    /// Materialization may use an _increased_ version of this if it truncates
    /// the alignment at the first stop codon encountered.
    pub rpad: usize,
}

impl Product<'_> {
    /// Returns the length of the full coding sequence as specified in the
    /// corresponding [`ProductSpec`].
    ///
    /// [`ProductSpec`]: crate::config::ProductSpec
    pub(crate) fn full_cds_len(&self) -> usize {
        self.lpad + self.cds_aligned_len + self.rpad
    }
}
