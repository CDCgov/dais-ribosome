use std::ops::Range;

/// Exon specification for a protein product (ctype stripped).
#[derive(Debug, Clone)]
pub(crate) struct Exons {
    /// The optionally required codon which must be present at the start of the
    /// alignment in order for the product to be included.
    ///
    /// This contains solely `ACGT`.
    pub(crate) required_start: Option<[u8; 3]>,

    /// The coordinates of the exons within the reference and coding sequence.
    ///
    /// The exons are ordered by `cds_range`, which form a partition of
    /// `0..cds_len` where `cds_len` is the total length of the coding sequence
    /// (a multiple of 3). The `ref_range` fields are in order, although up to 2
    /// nucleotides overlap is allowed between ranges. Note that any repeated
    /// indices are represented twice with distinct coordinates in the coding
    /// sequence.
    ///
    /// This vector is non-empty.
    pub(crate) coords: Vec<ExonCoords>,
}

impl Exons {
    /// The coordinates of the first exon.
    #[inline]
    #[allow(dead_code)]
    pub fn first(&self) -> &ExonCoords {
        self.coords.first().expect("The coords field of Exons should be non-empty")
    }

    /// The coordinates of the last exon.
    #[inline]
    pub fn last(&self) -> &ExonCoords {
        self.coords.last().expect("The coords field of Exons should be non-empty")
    }

    /// The length of the coding sequence as defined by the exons.
    ///
    /// This is guaranteed to be a multiple of 3.
    #[inline]
    pub fn cds_len(&self) -> usize {
        self.last().cds_range.end
    }
}

/// The coordinates of an exon (coding sequence) within a reference and coding
/// sequence.
///
/// ## Validity
///
/// The two ranges should be the same length and must be non-empty. They need
/// not be multiples of 3 in length.
#[derive(Debug, Clone)]
pub struct ExonCoords {
    /// The 0-based end-exclusive range where the exon occurs within the
    /// reference sequence.
    pub(crate) ref_range: Range<usize>,

    /// The 0-based end-exclusive range where the exon occurs within the coding
    /// sequence.
    pub(crate) cds_range: Range<usize>,
}
