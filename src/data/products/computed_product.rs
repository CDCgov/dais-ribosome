use crate::data::{products::ProductSpec, ranges::CdsStateRange};
use zoe::prelude::*;

/// Pre-computed product data ready for output.
#[derive(Debug)]
pub struct ComputedProduct<'a> {
    pub(crate) product_ranges: Vec<CdsStateRange>,
    pub(crate) product_spec:   &'a ProductSpec,

    /// CDS sequence (without deletions, includes insertions)
    pub cds_seq:         Nucleotides,
    /// CDS alignment (with `-` for deletions, no insertions)
    pub cds_aln:         Nucleotides,
    /// SHA1 hash of cleaned CDS sequence
    pub cds_id:          String,
    /// Amino acid sequence
    pub aa_seq:          AminoAcids,
    /// Amino acid alignment (with `-` for deletions)
    pub aa_aln:          AminoAcids,
    /// MD5 hash of cleaned AA sequence (variant hash)
    pub variant_hash:    String,
    /// Whether any insertion exists in this product
    pub has_insertion:   bool,
    /// Whether any insertion or deletion causes a frameshift (length % 3 != 0)
    pub has_shift_indel: bool,
    /// Query nucleotide coordinates (e.g., "1..45;48..90")
    pub query_coords:    String,
    /// CDS nucleotide coordinates (e.g., "1..45")
    pub cds_coords:      String,
    /// Computed insertions for this product
    pub insertions:      Vec<ComputedInsertion>,
    /// Computed deletions for this product
    pub deletions:       Vec<ComputedDeletion>,
}

impl ComputedProduct<'_> {
    pub(crate) fn trailing_cds_gap_len(&self) -> usize {
        // Use the furthest CDS end from any M or D state to avoid
        // double-padding when a deletion covers trailing positions.
        let last_cds_end = self
            .product_ranges
            .iter()
            .rev()
            .find_map(|s| match s {
                CdsStateRange::M(m) => Some(m.cds_range.end),
                CdsStateRange::D(d) => Some(d.cds_range.end),
                _ => None,
            })
            .unwrap_or(0);

        self.product_spec.exons.total_cds_length - last_cds_end
    }
}

/// Pre-computed insertion data ready for output.
#[derive(Debug)]
pub struct ComputedInsertion {
    /// Upstream amino acid position (1-based)
    pub upstream_aa:          usize,
    /// Upstream nucleotide position (1-based)
    pub upstream_nt:          usize,
    /// Inserted nucleotides
    pub inserted_nucleotides: Nucleotides,
    /// Inserted residues (translated)
    pub inserted_residues:    AminoAcids,
    /// Codon shift (0, 1, or 2)
    pub codon_shift:          usize,
    /// Whether this insertion should be filtered from AA output
    pub filtered:             bool,
}

impl ComputedInsertion {
    /// Create a `ComputedInsertion` from raw insertion data.
    ///
    /// `upstream_nt_pos` is the 1-based upstream nucleotide position in CDS space.
    pub fn new(upstream_nt_pos: usize, slice: &[u8]) -> Self {
        let ins_len = slice.len();
        let inserted_nucleotides = Nucleotides::from_vec_unchecked(slice.to_vec());

        let upstream_aa = upstream_nt_pos / 3;
        let codon_shift = upstream_nt_pos % 3;

        let (inserted_residues, filtered) = if ins_len < 3 {
            (AminoAcids::from_vec_unchecked(vec![b'?']), true)
        } else if slice.iter().all(|&b| b == b'n' || b == b'N') {
            (AminoAcids::from_vec_unchecked(vec![b'X']), true)
        } else {
            (inserted_nucleotides.translate_to_stop(), false)
        };

        ComputedInsertion {
            upstream_aa,
            upstream_nt: upstream_nt_pos,
            inserted_nucleotides,
            inserted_residues,
            codon_shift,
            filtered,
        }
    }
}

/// Pre-computed deletion data ready for output.
#[derive(Debug)]
pub struct ComputedDeletion {
    /// Deletion start in amino acid coordinates (1-based)
    pub del_aa_start:  usize,
    /// Deletion end in amino acid coordinates (1-based)
    pub del_aa_end:    usize,
    /// Deletion length in amino acids
    pub del_aa_len:    usize,
    /// Whether deletion is in-frame
    pub in_frame:      bool,
    /// Deletion start in CDS coordinates (1-based)
    pub del_cds_start: usize,
    /// Deletion end in CDS coordinates (1-based)
    pub del_cds_end:   usize,
    /// Deletion length in nucleotides
    pub del_cds_len:   usize,
}
