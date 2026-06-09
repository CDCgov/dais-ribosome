use crate::ranges::{InclusiveDisplay, InsertionIdx, RangeExt};
use std::{cmp::Ordering, ops::Range};

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
    /// Constructs a new [`Exons`] object from a vector of parsed coordinates.
    ///
    /// ## Validity
    ///
    /// If specified, `required_start` should contain solely characters in
    /// `ACGT`.
    ///
    /// ## Errors
    ///
    /// - The `coords` field must be non-empty.
    /// - The sum of all exon lengths must be a multiple of three.
    ///
    /// For each pair of consecutive exons:
    ///
    /// - The exons must be in the correct order, so that both endpoints of the
    ///   second exon are at least the corresponding endpoints of the first exon
    ///   (and at least one of the inequalities is strict).
    /// - The endpoint of the first exon cannot equal the starting point of the
    ///   next exon.
    /// - The overlap between the two exons can be at most 2.
    ///
    /// For each set of three consecutive exons:
    ///
    /// - The three exons cannot all share a single region of overlap.
    pub fn new(coords: Vec<ExonCoords>, required_start: Option<[u8; 3]>) -> std::io::Result<Self> {
        // Validate coords is non-empty and total length is multiple of 3
        let cds_len = coords
            .last()
            .ok_or(std::io::Error::other("The coordinates for the exon must be non-empty"))?
            .cds_range
            .end;

        if !cds_len.is_multiple_of(3) {
            return Err(std::io::Error::other(
                "The length of the coding sequence (sum of all exon lengths) was not a multiple of 3.",
            ));
        }

        let mut overlapped_regions = Vec::new();
        let mut noncoding_regions = Vec::new();

        for exons in coords.array_windows() {
            Self::validate_two_exons(exons)?;
            let [exon1, exon2] = exons;

            // Overlapping and containing a noncoding region between them are
            // mutually exclusive, so use "else if"
            if let Some(overlapped) = ExonOverlapCoords::new(exon1, exon2) {
                overlapped_regions.push(overlapped);
            } else if let Some(noncoding) = NoncodingCoords::new(exon1, exon2) {
                noncoding_regions.push(noncoding);
            }
        }

        for exons in coords.array_windows() {
            Self::validate_three_exons(exons)?;
        }

        let exons = Exons {
            required_start,
            coords,
            overlapped_regions,
            noncoding_regions,
        };

        Ok(exons)
    }

    /// A helper function for [`Exons::new`] that checks the requirements for
    /// two consecutive exons.
    ///
    /// ## Errors
    ///
    /// - The exons must be in the correct order, so that both endpoints of the
    ///   second exon are at least the corresponding endpoints of the first exon
    ///   (and at least one of the inequalities is strict).
    /// - The endpoint of the first exon cannot equal the starting point of the
    ///   next exon.
    /// - The overlap between the two exons can be at most 2.
    fn validate_two_exons(exons: &[ExonCoords; 2]) -> std::io::Result<()> {
        /// The maximum amount of overlap allowed between exons.
        ///
        /// SARS-CoV-2 requires -1 exon-to-exon frameshift with other viruses
        /// reported up to -2.
        const MAX_DUPLICATED_OVERLAP_NT: usize = 2;

        let [left, right] = exons;

        // Ensure correct ordering
        match left.ref_range.relaxed_cmp(&right.ref_range) {
            Some(Ordering::Less) => {}
            Some(Ordering::Greater) => {
                return Err(std::io::Error::other(format!(
                    "Exons out of order! Found {} then {}",
                    left.ref_range.display_inclusive(),
                    right.ref_range.display_inclusive(),
                )));
            }
            Some(Ordering::Equal) => {
                return Err(std::io::Error::other(format!(
                    "Found the same exon twice: {}",
                    left.ref_range.display_inclusive()
                )));
            }
            None => {
                return Err(std::io::Error::other(format!(
                    "One exon cannot completely contain another! Found {} then {}",
                    left.ref_range.display_inclusive(),
                    right.ref_range.display_inclusive(),
                )));
            }
        };

        // Prevent perfectly adjacent exons (there should either be overlap or
        // non-coding region)
        if left.ref_range.end == right.ref_range.start {
            return Err(std::io::Error::other(format!(
                "Two exons are perfectly adjacent, and should therefore be represented as a single exon. Found {} then {}",
                left.ref_range.display_inclusive(),
                right.ref_range.display_inclusive(),
            )));
        }

        // Prevent overlapping exons that overlap by more than
        // MAX_DUPLICATED_OVERLAP_NT
        //
        // Exclusive index - inclusive index is valid length
        let overlap_nt = left.ref_range.end.saturating_sub(right.ref_range.start);
        if overlap_nt > MAX_DUPLICATED_OVERLAP_NT {
            return Err(std::io::Error::other(format!(
                "Exon overlap exceeds {MAX_DUPLICATED_OVERLAP_NT} nt! Found {} then {}",
                left.ref_range.display_inclusive(),
                right.ref_range.display_inclusive(),
            )));
        }

        Ok(())
    }

    /// A helper function for [`Exons::new`] that checks the requirements for
    /// three consecutive exons.
    ///
    /// ## Errors
    ///
    /// The three exons cannot all share a single region of overlap.
    fn validate_three_exons(exons: &[ExonCoords; 3]) -> std::io::Result<()> {
        // Prevent a single region of overlap from involving more than 2
        // exons
        let [first, middle, last] = exons;

        let overlap_nt = first.ref_range.end.saturating_sub(last.ref_range.start);

        if overlap_nt > 0 {
            return Err(std::io::Error::other(format!(
                "A single region of overlap cannot involve more than 2 exons within a given protein product. Found {}, {}, then {}, which all overlap",
                first.ref_range.display_inclusive(),
                middle.ref_range.display_inclusive(),
                last.ref_range.display_inclusive(),
            )));
        }

        Ok(())
    }

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
