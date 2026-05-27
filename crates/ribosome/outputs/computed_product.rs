use crate::{
    data::ranges::InsertionIdx,
    ranges::{CdsCoord, CdsDeletionRange},
};
use std::ops::Range;
use zoe::prelude::*;

/// The three cases for materializing a product.
///
/// The normal output is [`MaybeComputedProduct::Ok`], which contains the
/// [`ComputedProduct`].
pub enum MaybeComputedProduct<'a> {
    /// The normal output for a non-empty product with at least one match state.
    Ok(ComputedProduct<'a>),
    /// An empty product, caused by the genome alignment not intersecting the
    /// exons.
    Empty(EmptyProduct<'a>),
    /// A fully-deleted product, caused by all of the exons being covered by
    /// deletions.
    Deleted(DeletedProduct<'a>),
}

/// An empty product, caused by the genome alignment not intersecting the exons.
#[derive(Debug)]
pub struct EmptyProduct<'a> {
    /// The protein product name (e.g., `HA`, `HA-signal`).
    pub name: &'a str,

    /// The total padding for `cds_aln`, which is equal to the length of the
    /// coding sequence from the specifications since the product is empty.
    pub cds_aln_pad: usize,
}

#[derive(Debug)]
pub struct DeletedProduct<'a> {
    /// The protein product name (e.g., `HA`, `HA-signal`).
    pub name:         &'a str,
    /// The aligned coding sequence for the protein (with `-` for deletions but
    /// no insertions).
    ///
    /// This will only contain uppercase IUPAC, padding `.`, and gaps `-`. Both
    /// `U` and `T` are allowed.
    pub cds_aln:      Nucleotides,
    /// The aligned amino acid sequence for the protein (with `-` for deletions
    /// but no insertions).
    ///
    /// This will only contain uppercase IUPAC, partial codons, stop codons,
    /// padding `.`, and gaps `-`.
    pub aa_aln:       AminoAcids,
    /// The computed deletions within the product.
    ///
    /// Unlike `insertions`, deletions of any length are included.
    pub deletion:     ComputedDeletion,
    /// The amount of right padding that occurs after `cds_aln`.
    ///
    /// Padding of this length can be added to the end of `cds_aln` to get an
    /// aligned sequence spanning the full coding sequence length. This does not
    /// include trailing deletions, since these are present in `cds_aln`
    /// already.
    pub cds_aln_rpad: usize,
}

/// A computed product, with materialized coding and amino acid sequences.
///
/// If an in-frame stop codon is encountered (and it does not appear in an
/// insertion), then all sequences in this struct will be truncated to exclude
/// any residues after it. Stop codons inside insertions do not truncate the
/// product, so if a stop codon appears in the middle of `aa_seq`, this is the
/// cause.
#[derive(Clone, Eq, PartialEq, Hash, Debug)]
pub struct ComputedProduct<'a> {
    /// The protein product name (e.g., `HA`, `HA-signal`).
    pub name:            &'a str,
    /// The unaligned coding sequence for the protein (with insertions but no
    /// deletions).
    ///
    /// This will only contain unaligned uppercase IUPAC. Both `U` and `T` are
    /// allowed.
    pub cds_seq:         Nucleotides,
    /// The aligned coding sequence for the protein (with `-` for deletions but
    /// no insertions).
    ///
    /// This will only contain uppercase IUPAC, left padding `.`, and gaps `-`.
    /// Both `U` and `T` are allowed. Right padding is not included, but the
    /// length of it can be obtained from [`ComputedProduct::cds_aln_rpad`].
    pub cds_aln:         Nucleotides,
    /// The SHA1 hash of the unaligned coding sequence (`cds_seq`).
    pub cds_id:          String,
    /// The unaligned amino acid sequence for the protein (with insertions but
    /// no deletions).
    ///
    /// This sequence is not a direct translation of `cds_seq`, which could
    /// produce noisy or unexpected results when a frameshift indel or an indel
    /// not aligning to codon boundaries is present. Instead, this is an effort
    /// to preserve the semantics of `aa_aln` but include insertions.
    /// Specifically, deletions from `aa_aln` are removed, and insertions are
    /// spliced into the amino acid sequence. Partial codons are introduced for
    /// frameshift indels (rather than altering the frame of translation), and
    /// indels in the middle of a codon are spliced before that codon's amino
    /// acid.
    ///
    /// This will only contain unaligned uppercase IUPAC, partial codons, and
    /// stop codons.
    pub aa_seq:          AminoAcids,
    /// The aligned amino acid sequence for the protein (with `-` for deletions
    /// but no insertions).
    ///
    /// This will only contain uppercase IUPAC, partial codons, stop codons,
    /// left padding `.`, and gaps `-`. Right padding is not included, but the
    /// length of it can be obtained from [`ComputedProduct::aa_aln_rpad`].
    pub aa_aln:          AminoAcids,
    /// The MD5 hash of the unaligned amino acid sequence (`aa_seq`).
    pub variant_hash:    String,
    /// Whether any insertion exists in this product.
    ///
    /// This includes insertions that were filtered, so it is possible that
    /// `has_insertions` is true but `insertions` is empty. A stop extension (if
    /// present) is also counted as an insertion.
    pub has_insertion:   bool,
    /// Whether any insertion or deletion causes a frameshift (i.e., the length
    /// is not divisible by 3).
    ///
    /// Flanking indels at either end of the product are not counted.
    pub has_shift_indel: bool,
    /// The coordinates within the original query that were used to form the
    /// unaligned `cds_seq`.
    ///
    /// After sanitizing the original query, slicing the original query at these
    /// indices and concatenating the results will yield `cds_seq`.
    pub query_coords:    Vec<Range<usize>>,
    /// The coding sequence coordinates corresponding to the
    /// `query_coordinates`.
    ///
    /// Each coding sequence coordinate is either a range (corresponding to a
    /// matched region) or an [`InsertionIdx`] (corresponding to an insertion in
    /// the query). Discontinuities in the ranges imply a deletion in the query.
    /// This vector will be the same length as `query_coords`.
    pub cds_coords:      Vec<CdsCoord>,
    /// The computed insertions within the product.
    ///
    /// Only insertions that have a length at least 3 are included. Insertions
    /// solely containing `N` are excluded. If it meets the criteria, a stop
    /// extension is included as an insertion.
    pub insertions:      Vec<ComputedInsertion>,
    /// The computed deletions within the product.
    ///
    /// Unlike `insertions`, deletions of any length are included.
    pub deletions:       Vec<ComputedDeletion>,
    /// The amount of right padding that occurs after `cds_aln`.
    ///
    /// Padding of this length can be added to the end of `cds_aln` to get an
    /// aligned sequence spanning the full coding sequence length. This does not
    /// include trailing deletions, since these are present in `cds_aln`
    /// already.
    pub cds_aln_rpad:    usize,
}

impl DeletedProduct<'_> {
    /// The amount of right padding that occurs after `aa_aln`.
    ///
    /// Padding of this length can be added to the end of `aa_aln` to get an
    /// aligned sequence spanning the full coding protein length. This does not
    /// include trailing deletions, since these are present in `aa_aln` already.
    pub fn aa_aln_rpad(&self) -> usize {
        // Floor divide since any leftover bases are already accounted for with
        // a partial codon
        self.cds_aln_rpad / 3
    }
}

impl ComputedProduct<'_> {
    /// The amount of right padding that occurs after `aa_aln`.
    ///
    /// Padding of this length can be added to the end of `aa_aln` to get an
    /// aligned sequence spanning the full coding protein length. This does not
    /// include trailing deletions, since these are present in `aa_aln` already.
    pub fn aa_aln_rpad(&self) -> usize {
        // Floor divide since any leftover bases are already accounted for with
        // a partial codon
        self.cds_aln_rpad / 3
    }
}

/// A computed insertion, with materialized nucleotide and amino acid sequences.
///
/// Even if a stop codon occurs within the insertion, the sequence fields will
/// contain the full insertion.
#[derive(Clone, Eq, PartialEq, Hash, Debug)]
pub struct ComputedInsertion {
    /// The upstream amino acid position (1-based), which is the position
    /// _after_ which the insertion occurs.
    ///
    /// If the insertion interrupts a codon (i.e., `codon_shift` is nonzero),
    /// then this is rounded down (the insertion is treated as occuring before
    /// the split codon). This means that this field may be 0, which would
    /// represent an insertion within the first codon.
    pub upstream_aa_pos: usize,
    /// The upstream nucleotide position (1-based), which is the position
    /// _after_ which the insertion occurs.
    pub upstream_nt_pos: usize,
    /// The inserted nucleotides.
    ///
    /// This will only contain unaligned uppercase IUPAC. Both `U` and `T` are
    /// allowed.
    pub inserted_nt:     Nucleotides,
    /// A direct translation of `inserted_nt` to amino acids.
    ///
    /// This will only contain unaligned uppercase IUPAC, partial codons, and
    /// stop codons. A partial codon `~` is added to the end if the length is
    /// not a multiple of 3.
    pub inserted_aa:     AminoAcids,
    /// The codon shift of the insertion, which is the number of bases between
    /// the last codon and the insertion.
    ///
    /// 0 means that the insertion occurs between codons. 1 means that the
    /// insertion appears after the first base of a codon. 2 means that the
    /// insertion appears after the second base of a codon.
    pub codon_shift:     usize,
}

impl ComputedInsertion {
    /// Creates a [`ComputedInsertion`] from raw insertion data, and returns
    /// whether it should be filtered.
    ///
    /// The `cds_index` argument is the 0-based upstream nucleotide index in the
    /// coding sequence where the insertion occurs. The second return argument
    /// is true (indicating that it should be filtered) if the insertion length
    /// is less than 3 or the insertion is all `N`.
    ///
    /// ## Validity
    ///
    /// The `slice` must contain unaligned, uppercase IUPAC. It may contain `U`
    /// or `T`. This is true for any slice of [`QueryRecord::nucleotides`].
    ///
    /// [`QueryRecord::nucleotides`]: crate::data::QueryRecord::nucleotides
    pub(crate) fn new(cds_index: InsertionIdx, slice: &[u8]) -> (Self, bool) {
        let ins_len = slice.len();
        let inserted_nt = Nucleotides::from(slice);

        let aa_insertion_idx = cds_index.to_aa_idx();
        let codon_shift = cds_index.codon_shift();

        let (inserted_aa, filtered) = if ins_len < 3 || slice.iter().all(|&b| b == b'N') {
            // Do not include the all N insertion or shorter than 3 insertions
            // in the unaligned amino acid sequence output
            (AminoAcids::new(), true)
        } else {
            (inserted_nt.translate(), false)
        };

        (
            ComputedInsertion {
                upstream_aa_pos: aa_insertion_idx.left_pos(),
                upstream_nt_pos: cds_index.left_pos(),
                inserted_nt,
                inserted_aa,
                codon_shift,
            },
            filtered,
        )
    }
}

/// A computed deletion.
#[derive(Clone, Eq, PartialEq, Hash, Debug)]
pub struct ComputedDeletion {
    /// The start position of the deletion in amino acid coordinates (1-based,
    /// inclusive).
    pub del_aa_start:  usize,
    /// The end position of the deletion in amino acid coordinates (1-based,
    /// inclusive).
    pub del_aa_end:    usize,
    /// The deletion length in amino acids.
    ///
    /// This is equal to `del_cds_len.div_ceil(3)`, which can be interpreted as
    /// the number of codons-worth of nucleotides deleted. This is not
    /// necessarily equal to the number of codons intersected by the deletion,
    /// such as a two-base deletion spanning the end of one codon and the
    /// beginning of another.
    pub del_aa_len:    usize,
    /// Whether deletion is in-frame (both the CDS start position and length
    /// must be multiples of 3).
    pub in_frame:      bool,
    /// The start position of the deletion in coding sequence coordinates
    /// (1-based, inclusive).
    pub del_cds_start: usize,
    /// The end position of the deletion in coding sequence coordinates
    /// (1-based, inclusive).
    pub del_cds_end:   usize,
    /// The deletion length in nucleotides.
    pub del_cds_len:   usize,
}

impl ComputedDeletion {
    /// Creates a [`ComputedDeletion`] from raw deletion data.
    pub(crate) fn new(del: &CdsDeletionRange) -> Self {
        let in_frame = del.cds_range.start.is_multiple_of(3) && del.len().is_multiple_of(3);

        // The 1-based inclusive start. Floor divide so that
        // deleting any part of a codon deletes the amino acid. Add
        // 1 to make it 1-based.
        let del_aa_start = (del.cds_range.start / 3) + 1;

        // The 1-based inclusive end. Ceiling divide so that
        // deleting any part of a codon deletes the amino acid.
        // 0-based exclusive end is equivalent to 1-based inclusive
        // end.
        let del_aa_end = del.cds_range.end.div_ceil(3);

        let del_cds_len = del.len();
        let del_aa_len = del_cds_len.div_ceil(3);

        ComputedDeletion {
            del_aa_start,
            del_aa_end,
            del_aa_len,
            in_frame,
            del_cds_start: del.cds_range.start + 1,
            del_cds_end: del.cds_range.end,
            del_cds_len,
        }
    }
}
