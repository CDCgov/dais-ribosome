//! Key types used throughout the crate for indexing maps.

use std::fmt;

/// Key for reference sequences: combines reference ID and compound type.
///
/// Used to index into `ReferenceMap` to retrieve reference sequences.
#[derive(Hash, PartialEq, Eq, Debug, Clone)]
pub struct RefKey {
    pub reference_id:  String,
    pub compound_type: String,
}

impl RefKey {
    /// Create a new RefKey.
    pub fn new(reference_id: impl Into<String>, compound_type: impl Into<String>) -> Self {
        Self {
            reference_id:  reference_id.into(),
            compound_type: compound_type.into(),
        }
    }

    /// Parse a pipe-delimited name: `reference_id|compound_type`
    pub fn parse(name: &str) -> Option<Self> {
        let mut parts = name.split('|');
        let reference_id = parts.next()?;
        let compound_type = parts.next()?;
        Some(Self::new(reference_id, compound_type))
    }
}

impl From<(String, String)> for RefKey {
    fn from((ref_id, ctype): (String, String)) -> Self {
        Self {
            reference_id:  ref_id,
            compound_type: ctype,
        }
    }
}

impl fmt::Display for RefKey {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}|{}", self.reference_id, self.compound_type)
    }
}

/// Key for CDS specifications: combines reference ID and protein product name.
///
/// Used to index into `CdsSpecMap` to retrieve exon coordinates.
#[derive(Hash, PartialEq, Eq, Debug, Clone)]
pub struct SpecKey {
    pub reference_id:    String,
    pub protein_product: String,
}

impl SpecKey {
    /// Create a new SpecKey.
    pub fn new(reference_id: impl Into<String>, protein_product: impl Into<String>) -> Self {
        Self {
            reference_id:    reference_id.into(),
            protein_product: protein_product.into(),
        }
    }
}

impl From<(&str, &str)> for SpecKey {
    fn from((ref_id, prot): (&str, &str)) -> Self {
        SpecKey::new(ref_id, prot)
    }
}

impl fmt::Display for SpecKey {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}|{}", self.reference_id, self.protein_product)
    }
}

/// Key for codon position weights: combines position and codon triplet.
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub struct CodonKey {
    pub position: u32,
    pub codon:    [u8; 3],
}

impl CodonKey {
    /// Create a new CodonKey, normalizing the codon to uppercase.
    pub fn new(position: u32, codon: [u8; 3]) -> Self {
        let mut key = Self { position, codon };
        key.codon.make_ascii_uppercase();
        key
    }
}

impl fmt::Debug for CodonKey {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}\t{}{}{}",
            self.position, self.codon[0] as char, self.codon[1] as char, self.codon[2] as char,
        )
    }
}

impl fmt::Display for CodonKey {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}:{}{}{}",
            self.position, self.codon[0] as char, self.codon[1] as char, self.codon[2] as char,
        )
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

    #[test]
    fn test_codon_key_uppercase() {
        let key = CodonKey::new(42, *b"atg");
        assert_eq!(key.codon, *b"ATG");
    }
}
