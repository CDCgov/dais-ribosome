use std::ops::Range;

/// Exon specification for a protein product (ctype stripped).
#[derive(Debug, Clone)]
pub(crate) struct Exons {
    pub(crate) required_start:   Option<[u8; 3]>,
    /// The coordinates and offsets of the exons (coding sequences) within the
    /// reference.
    ///
    /// These are guaranteed to be in order and non-overlapping.
    pub(crate) coords:           Vec<ExonCoords>,
    /// The total length of all the exons for this entry in
    /// the `cds-spec.tsv` file.
    pub(crate) total_cds_length: usize,
}

/// The coordinates of an exon (coding sequence) within a reference, as well as
/// an offset from the previous exon.
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

    /// The offset of the reference coordinates to the coding sequence
    /// coordinates (i.e., the number of non-coding residues up until this
    /// exon).
    ///
    /// This is subtracted from reference coordinates to get coding sequence
    /// coordinates.
    pub(crate) ref_to_cds_offset: usize,
}
