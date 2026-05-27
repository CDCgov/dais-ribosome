//! Sequence hashing utilities ported from `Omics::Nucleotide` and
//! `Omics::AminoAcids` in `ifx-perl`

use zoe::data::{nucleotides::Nucleotides, types::amino_acids::AminoAcids};

/// Computes the SHA1 hex of the unaligned nucleotide sequence.
///
/// This is a wrapper around [`sha1_hex`].
///
/// ## Validity
///
/// This must only contain unaligned uppercase IUPAC. Both `U` and `T` are
/// allowed. `seq` must be non-empty.
#[inline]
pub(crate) fn nt_id_iupac(seq: &Nucleotides) -> String {
    sha1_hex(seq.as_bytes())
}

/// Returns MD5 hex of the unaligned amino acid sequence.
///
/// This is a wrapper around [`md5_hex`].
///
/// ## Validity
///
/// This must only contain unaligned uppercase IUPAC, partial codons, and stop
/// codons. `seq` must be non-empty.
#[inline]
pub(crate) fn variant_hash_iupac(seq: &AminoAcids) -> String {
    md5_hex(seq.as_bytes())
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
    fn test_variant_hash_basic() {
        let seq: AminoAcids = b"ACDE".into();
        let hash = variant_hash_iupac(&seq);

        assert_eq!(hash.len(), 32);
        assert!(hash.chars().all(|c| c.is_ascii_hexdigit()));
    }

    #[test]
    fn test_variant_hash_preserves_tilde() {
        let with_tilde: AminoAcids = b"ACDE~".into();
        let without_tilde: AminoAcids = b"ACDE".into();

        assert_ne!(variant_hash_iupac(&with_tilde), variant_hash_iupac(&without_tilde));
    }
}
