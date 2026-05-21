//! Key types used throughout the crate for indexing maps.

use std::fmt;

/// A key for reference sequences, combining a reference ID and a compound type.
#[derive(Hash, PartialEq, Eq, Debug, Clone)]
pub struct RefKey {
    /// The reference ID of the reference sequence (e.g., `ANHUI01`,
    /// `PHUKET3073`).
    pub reference_id:  String,
    /// The compound type of the reference sequence (e.g., `A_HA_H7`, `B_HA`).
    pub compound_type: String,
}

impl RefKey {
    /// Creates a new [`RefKey`] by combining a `reference_id` (e.g., `ANHUI01`)
    /// and `compound_type` (e.g., `A_HA_H7`).
    pub fn new(reference_id: impl Into<String>, compound_type: impl Into<String>) -> Self {
        Self {
            reference_id:  reference_id.into(),
            compound_type: compound_type.into(),
        }
    }

    /// Parses a pipe-delimited FASTA header of the form
    /// `reference_id|compound_type`.
    ///
    /// Any additional pipe-delimited fields are ignored.
    pub fn parse(name: &str) -> Option<Self> {
        let mut parts = name.split('|');
        let reference_id = parts.next()?;
        let compound_type = parts.next()?;
        Some(Self::new(reference_id, compound_type))
    }
}

impl fmt::Display for RefKey {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}|{}", self.reference_id, self.compound_type)
    }
}

/// A key for coding sequence (CDS) specifications, combining a reference ID and
/// protein product name.
///
/// Used to index into [`CdsSpecMap`] to retrieve exon coordinates.
///
/// [`CdsSpecMap`]: crate::config::cds_spec::CdsSpecMap
#[derive(Clone, Eq, PartialEq, Hash, Debug)]
pub struct SpecKey {
    /// The reference ID of the reference sequence (e.g., `A_HA_H7`, `B_HA`).
    pub reference_id: String,
    /// The protein product name (e.g., `HA`, `HA-signal`).
    pub product_name: String,
}

impl SpecKey {
    /// Creates a new [`SpecKey`] by combining a `reference_id` (e.g.,
    /// `ANHUI01`) and `protein_product` (e.g., `HA`).
    pub fn new(reference_id: impl Into<String>, product_name: impl Into<String>) -> Self {
        Self {
            reference_id: reference_id.into(),
            product_name: product_name.into(),
        }
    }
}

impl fmt::Display for SpecKey {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}|{}", self.reference_id, self.product_name)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ref_key_parse() {
        let key = RefKey::parse("A_HA|H3N2").unwrap();
        assert_eq!(key.reference_id, "A_HA");
        assert_eq!(key.compound_type, "H3N2");
    }

    #[test]
    fn test_ref_key_parse_missing_pipe() {
        assert!(RefKey::parse("no_pipe").is_none());
    }
}
