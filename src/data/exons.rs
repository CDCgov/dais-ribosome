use std::ops::Range;

/// Exon specification for a protein product (ctype stripped).
#[derive(Debug, Clone)]
pub(crate) struct Exons {
    pub(crate) required_start:   Option<[u8; 3]>,
    /// The coordinates of the exons (coding sequences) within the reference
    /// and coding sequence.
    ///
    /// These are guaranteed to be in order. Reference overlaps of up to two
    /// nucleotides are allowed and are duplicated in CDS order.
    pub(crate) coords:           Vec<ExonCoords>,
    /// The total length of all the exons for this entry in
    /// the `cds-spec.tsv` file.
    pub(crate) total_cds_length: usize,
}

/// The coordinates of an exon (coding sequence) within a reference and coding
/// sequence.
#[derive(Debug, Clone)]
pub struct ExonCoords {
    /// The 0-based end-exclusive range where the exon occurs within the
    /// reference sequence.
    ///
    /// ## Validity
    ///
    /// This range is always non-empty. This range need not be a multiple of 3
    /// in length.
    pub(crate) ref_range: Range<usize>,

    /// The 0-based end-exclusive range where the exon occurs within the coding
    /// sequence.
    ///
    /// This range always has the same length as [`Self::ref_range`].
    pub(crate) cds_range: Range<usize>,
}
