use std::ops::Range;

/// Exon specification for a protein product (ctype stripped).
#[derive(Debug, Clone)]
pub(crate) struct Exons {
    pub(crate) required_start:   Option<[u8; 3]>,
    /// The coordinates and offsets of the exons (coding sequences) within the
    /// reference.
    pub(crate) coords:           Vec<ExonCoords>,
    /// The total length of all the exons for this entry in
    /// the `cds-spec.tsv` file.
    pub(crate) total_cds_length: usize,
}

/// The coordinates of an exon (coding sequence) within a reference, as well as
/// an offset from the previous exon.
#[derive(Debug, Clone)]
pub struct ExonCoords {
    /// The range where the exon occurs within the reference sequence.
    ///
    /// ## Validity
    ///
    /// This range is always non-empty.
    pub(crate) ref_range:         Range<usize>,
    /// The offset of this exon from the previous exon (i.e., the length of the
    /// intron between them).
    pub(crate) ref_to_cds_offset: usize,
}
