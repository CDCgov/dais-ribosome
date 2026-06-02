//! Sequence hashing utilities ported from `Omics::Nucleotide` and
//! `Omics::AminoAcids` in `ifx-perl`

use zoe::{
    data::{nucleotides::Nucleotides, types::amino_acids::AminoAcids},
    prelude::Len,
};

/// Computes the SHA1 hex of the unaligned nucleotide sequence.
///
/// ## Validity
///
/// This must only contain unaligned uppercase IUPAC. Both `U` and `T` are
/// allowed.
#[inline]
pub(crate) fn nt_id_iupac(seq: &Nucleotides) -> Option<String> {
    (!seq.is_empty()).then(|| sha1_hex(seq.as_bytes()))
}

/// Returns MD5 hex of uppercased sequence with `\n\r\t :.-` removed (keeps
/// `~`).
///
/// ## Validity
///
/// This must only contain unaligned uppercase IUPAC, partial codons, and stop
/// codons.
#[inline]
pub(crate) fn variant_hash_iupac(seq: &AminoAcids) -> Option<String> {
    (!seq.is_empty()).then(|| md5_hex(seq.as_bytes()))
}

/// Compute SHA1 hex digest of byte slice.
///
/// The string will always be a 40 character hexadecimal number.
#[inline]
fn sha1_hex(data: &[u8]) -> String {
    sha1_smol::Sha1::from(data).digest().to_string()
}

/// Compute MD5 hex digest of byte slice.
///
/// The string will always be a 32 character hexadecimal number.
#[inline]
fn md5_hex(data: &[u8]) -> String {
    // By default, [`md5`] automatically pads with zeros when using LowerHex
    // formatting.
    format!("{:x}", md5::compute(data))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_nt_id_empty() {
        let seq = Nucleotides::new();
        assert!(nt_id_iupac(&seq).is_none());
    }

    #[test]
    fn test_variant_hash_basic() {
        let seq: AminoAcids = b"ACDE".into();
        let hash = variant_hash_iupac(&seq).unwrap();

        assert_eq!(hash.len(), 32);
        assert!(hash.chars().all(|c| c.is_ascii_hexdigit()));
    }

    #[test]
    fn test_variant_hash_empty() {
        let seq = AminoAcids::new();
        assert!(variant_hash_iupac(&seq).is_none());
    }

    #[test]
    fn test_variant_hash_preserves_tilde() {
        let with_tilde: AminoAcids = b"ACDE~".into();
        let without_tilde: AminoAcids = b"ACDE".into();

        assert_ne!(variant_hash_iupac(&with_tilde), variant_hash_iupac(&without_tilde));
    }
}
