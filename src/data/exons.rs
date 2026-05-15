use std::ops::Range;

use crate::ranges::{InsertionIdx, RangeExt};

/// Exon specification for a protein product (ctype stripped).
#[derive(Clone, Debug)]
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

    /// The coordinates where any exons overlap, in a precomputed list to aid in
    /// indel shifting code.
    pub(crate) overlapped_regions: Vec<ExonOverlapCoords>,

    /// The coordinates of any non-coding regions between exons, in a
    /// precomputed list to aid in indel shifting code.
    pub(crate) noncoding_regions: Vec<NoncodingCoords>,
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

/// The coordinates of an overlap between two exons within a reference and
/// coding sequence.
///
/// ## Validity
///
/// The ranges must be the same length and must be non-empty. The ranges can be
/// of at most length 2 (the maximum allowed overlap).
#[derive(Debug, Clone)]
pub struct ExonOverlapCoords {
    /// The 0-based end-exclusive range where the overlap occurs within the
    /// reference sequence.
    #[allow(dead_code)]
    pub(crate) ref_range: Range<usize>,

    /// The 0-based end-exclusive range where the overlap occurs within the
    /// first exon.
    pub(crate) cds_range1: Range<usize>,

    /// The 0-based end-exclusive range where the overlap occurs within the
    /// second exon.
    pub(crate) cds_range2: Range<usize>,
}

impl ExonOverlapCoords {
    /// Identifies whether there is overlap between consecutive exons `exon1`
    /// and `exon2`, returning the overlap if it is present.
    pub fn new(exon1: &ExonCoords, exon2: &ExonCoords) -> Option<Self> {
        if let Some(ref_range) = exon1.ref_range.intersect_ordered(&exon2.ref_range) {
            let overlap_len = ref_range.len();
            let cds_range1 = exon1.cds_range.end - overlap_len..exon1.cds_range.end;
            let cds_range2 = exon2.cds_range.start..exon2.cds_range.start + overlap_len;

            Some(Self {
                ref_range,
                cds_range1,
                cds_range2,
            })
        } else {
            None
        }
    }

    pub fn cds_range(&self) -> Range<usize> {
        self.cds_range1.start..self.cds_range2.end
    }
}

/// The coordinates of a non-coding region between two exons within a reference
/// and coding sequence.
///
/// ## Validity
///
/// The range must be non-empty.
#[derive(Debug, Clone)]
pub struct NoncodingCoords {
    /// The 0-based end-exclusive range where the overlap occurs within the
    /// reference sequence.
    #[allow(dead_code)]
    pub(crate) ref_range: Range<usize>,

    /// The index of the noncoding region within the coding sequence. This can
    /// be viewed similarly to an insertion, hence reuse of [`InsertionIdx`].
    pub(crate) cds_index: InsertionIdx,
}

impl NoncodingCoords {
    /// Identifies whether there is a non-coding region between consecutive
    /// exons `exon1` and `exon2`, returning the coordinates if it is present.
    pub fn new(exon1: &ExonCoords, exon2: &ExonCoords) -> Option<Self> {
        debug_assert_eq!(exon1.cds_range.end, exon2.cds_range.start);

        if exon2.ref_range.start > exon1.ref_range.end {
            Some(Self {
                ref_range: exon1.ref_range.end..exon2.ref_range.start,
                cds_index: InsertionIdx::from_right_idx(exon2.cds_range.start),
            })
        } else {
            None
        }
    }
}
