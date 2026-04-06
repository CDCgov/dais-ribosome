use crate::{
    data::{
        ComputedProduct,
        products::{ProductSpec, incremental_products::ComputedIncrementalProducts},
        ranges::CdsStateRange,
    },
    hashing::{nt_id, variant_hash},
};
use std::ops::Range;
use zoe::prelude::*;

/// The aligned ranges for a single query against a single reference, using the
/// exons for one of the protein products.
///
/// Many [`Product`] values are stored in [`GenomeAndProductStates`], which
/// holds the alignments for all the protein products. Many of these are stored
/// in [`RibosomeOutput`], which holds the alignments for all reference IDs.
/// Many [`RibosomeOutput`] are generated and written in `main.rs` for each
/// query.
///
/// [`GenomeAndProductStates`]: crate::data::outputs::GenomeAndProductStates
/// [`RibosomeOutput`]: crate::data::outputs::RibosomeOutput
#[derive(Debug)]
pub(crate) struct Product<'a> {
    /// The information for the protein product being aligned against, including
    /// the name and exons.
    pub(crate) product_spec: &'a ProductSpec,

    /// The ranges within the exons that the query covers. This is initially
    /// formed by intersecting the query ranges with the exon ranges, then is
    /// tweaked.
    ///
    /// This is guaranteed to contain ordered and non-overlapping ranges. It
    /// does not begin or end with [`CdsStateRange::I`]. Within the aligned
    /// portion of the coding sequence (the exons), the ranges will be adjacent
    /// (forming a partition). However, there may be exons at the beginning or
    /// end which are not aligned against (or partially aligned against).
    pub(crate) product_ranges: Vec<CdsStateRange>,

    /// The number of bases in the coding sequence that were not aligned against
    /// in the beginning (i.e., not included in `product_ranges`).
    pub(crate) leading_cds_unaligned: usize,

    /// The number of bases in the coding sequence that were not aligned against
    /// at the end (i.e., not included in `product_ranges`).
    ///
    /// Materialization may use an _increased_ version of this if it truncates
    /// the alignment at the first stop codon encountered.
    pub(crate) trailing_cds_unaligned: usize,

    /// If this product's last exon ends at the stop extension position, this
    /// holds the query range of the stop extension nucleotides.
    pub(crate) stop_extension_query_range: Option<Range<usize>>,
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
                let num_to_consume = insertion.upstream_aa_pos - num_consumed;

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
