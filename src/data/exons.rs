use std::ops::Range;

/// Exon specification with a compound type, used during loading.
///
/// This contains the fields exactly as they were parsed from the
/// `cds-specs.tsv` file.
#[derive(Debug, Clone)]
pub(crate) struct CtypeExons {
    pub(crate) ctype:          String,
    pub(crate) required_start: Option<[u8; 3]>,
    pub(crate) coords:         Vec<ExonCoords>,
}

impl CtypeExons {
    /// Consumes `self`, splitting it into the compound type and [`Exons`].
    pub fn into_ctype_exons(self) -> (String, Exons) {
        let total_cds_length: usize = self.coords.iter().map(|r| r.ref_range.len()).sum();

        // TODO: This is not guaranteed to be true, so this should be an actual
        // check and result in an actual error. Does each range need to be a
        // multiple of three (e.g., parse_coordinate_ranges does validation)? Or
        // should the check be done here?
        debug_assert!(
            total_cds_length.is_multiple_of(3),
            "{} product was not in-frame: {total_cds_length}",
            self.ctype
        );

        (
            self.ctype,
            Exons {
                required_start: self.required_start,
                coords: self.coords,
                total_cds_length,
            },
        )
    }
}

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
