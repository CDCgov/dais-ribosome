use std::ops::Range;

use crate::{data::ranges::InsertionIdx, ranges::CdsCoord};
use zoe::prelude::*;

/// Pre-computed product data ready for output.
#[derive(Debug)]
pub struct ComputedProduct<'a> {
    /// The product name
    pub product_name:           &'a str,
    /// CDS sequence (without deletions, includes insertions)
    pub cds_seq:                Nucleotides,
    /// CDS alignment (with `-` for deletions, no insertions)
    pub cds_aln:                Nucleotides,
    /// SHA1 hash of cleaned coding sequence, or `None` if no DNA data remained
    /// after filtering.
    pub cds_id:                 Option<String>,
    /// Amino acid sequence
    pub aa_seq:                 AminoAcids,
    /// Amino acid alignment (with `-` for deletions)
    pub aa_aln:                 AminoAcids,
    /// MD5 hash of the cleaned amino acid sequence (variant hash), or `None` if
    /// no amino acid data remained after filtering.
    pub variant_hash:           Option<String>,
    /// Whether any insertion exists in this product
    pub has_insertion:          bool,
    /// Whether any insertion or deletion causes a frameshift (length % 3 != 0)
    pub has_shift_indel:        bool,
    /// Query nucleotide coordinates
    pub query_coords:           Vec<Range<usize>>,
    /// CDS nucleotide coordinates
    pub cds_coords:             Vec<CdsCoord>,
    /// Computed non-filtered insertions for this product
    pub insertions:             Vec<ComputedInsertion>,
    /// Computed deletions for this product
    pub deletions:              Vec<ComputedDeletion>,
    /// The number of unaligned bases at the end of the coding sequence that
    /// were soft clipped or appeared after the first stop codon.
    ///
    /// This does not include trailing deletions, so that this field can be used
    /// to render right padding without double counting deletions.
    pub trailing_cds_unaligned: usize,
}

/// Pre-computed insertion data ready for output.
#[derive(Debug)]
pub struct ComputedInsertion {
    /// The upstream amino acid position (1-based) after which the insertion
    /// occurs.
    ///
    /// If the insertion induces a frame shift (i.e., `codon_shift` is nonzero),
    /// then this is rounded down (the insertion is treated as occuring before
    /// the split codon). This means that this field may be 0, which would
    /// represent an insertion within the first codon (before the first amino
    /// acid).
    pub upstream_aa_pos:      usize,
    /// The upstream nucleotide position (1-based) after which the insertion
    /// occurs.
    pub upstream_nt_pos:      usize,
    /// All the inserted nucleotides, including any after a stop codon in
    /// `inserted_residues`.
    pub inserted_nucleotides: Nucleotides,
    /// Inserted residues, translated to [`AminoAcids`].
    ///
    /// This will only contain data up until the first stop codon, and may
    /// contain a partial codon if the insertion is of a length not divisible by
    /// 3.
    pub inserted_residues:    AminoAcids,
    /// The codon shift (0, 1, or 2).
    ///
    /// If nonzero, this indicates that the insertion causes a frameshift
    /// mutation.
    pub codon_shift:          usize,
}

impl ComputedInsertion {
    /// Creates a [`ComputedInsertion`] from raw insertion data, and returns
    /// whether it should be filtered.
    ///
    /// The `cds_index` argument is the 0-based upstream nucleotide position in
    /// the coding sequence where the insertion occurs. The second return
    /// argument is true (indicating that it should be filtered) if the
    /// insertion length is less than 3 or the insertion is all `N`.
    ///
    /// ## Validity
    ///
    /// The slice of the query range representing the insertion should contain
    /// unaligned, uppercase IUPAC bases.
    pub(crate) fn new(cds_index: InsertionIdx, slice: &[u8]) -> (Self, bool) {
        let ins_len = slice.len();
        let inserted_nucleotides = Nucleotides::from(slice);

        let aa_insertion_idx = cds_index.to_aa_idx();
        let codon_shift = cds_index.codon_shift();

        let (inserted_residues, filtered) = if ins_len < 3 || slice.iter().all(|&b| b == b'N') {
            // Do not include the all N insertion or shorter than 3 insertions
            // in the unaligned amino acid sequence output
            (AminoAcids::from(Vec::new()), true)
        } else {
            (inserted_nucleotides.translate(), false)
        };

        (
            ComputedInsertion {
                upstream_aa_pos: aa_insertion_idx.left_pos(),
                upstream_nt_pos: cds_index.left_pos(),
                inserted_nucleotides,
                inserted_residues,
                codon_shift,
            },
            filtered,
        )
    }
}

/// Pre-computed deletion data ready for output.
#[derive(Debug)]
pub struct ComputedDeletion {
    /// The start position of the deletion in amino acid coordinates (1-based,
    /// inclusive).
    pub del_aa_start:  usize,
    /// The end position of the deletion in amino acid coordinates (1-based,
    /// inclusive).
    pub del_aa_end:    usize,
    /// The deletion length in amino acids.
    pub del_aa_len:    usize,
    /// Whether deletion is in-frame (both the start position and length must be
    /// multiples of 3).
    pub in_frame:      bool,
    /// Deletion start in CDS coordinates (1-based).
    pub del_cds_start: usize,
    /// Deletion end in CDS coordinates (1-based).
    pub del_cds_end:   usize,
    /// Deletion length in nucleotides.
    pub del_cds_len:   usize,
}
