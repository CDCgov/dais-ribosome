use std::ops::Range;

/// Exon specification with ctype, used during loading.
#[derive(Debug, Clone)]
pub(crate) struct CtypeExons {
    pub(crate) ctype:          String,
    pub(crate) required_start: Option<[u8; 3]>,
    pub(crate) coords:         Vec<ExonCoords>,
}

impl CtypeExons {
    /// Consume self and split into (ctype, Exons).
    pub fn into_ctype_exons(self) -> (String, Exons) {
        let total_cds_length: usize = self.coords.iter().map(|r| r.ref_range.len()).sum();

        debug_assert!(
            total_cds_length.is_multiple_of(3),
            "{} product was not in-frame: {total_cds_length}",
            &self.ctype
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
    pub(crate) coords:           Vec<ExonCoords>,
    pub(crate) total_cds_length: usize,
}

#[derive(Debug, Clone)]
pub struct ExonCoords {
    pub(crate) ref_range:         Range<usize>,
    pub(crate) ref_to_cds_offset: usize,
}
