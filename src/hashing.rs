//! Sequence hashing utilities ported from `Omics::Nucleotide` and
//! `Omics::AminoAcids` in `ifx-perl`

use zoe::data::{ByteMap, RetainSequence, nucleotides::Nucleotides, types::amino_acids::AminoAcids};

/// A transformation mapping for use with [`retain_by_recoding`] which
/// uppercases amino acids and filters `\n\r\t :.-` (keeps `~` for partial
/// codons).
///
/// [`retain_by_recoding`]: RetainSequence::retain_by_recoding
const AA_HASH_FILTER: ByteMap = ByteMap::identity()
    .map_to_one(b"\n\r\t :.-", 0)
    .map_range(b'a'..=b'z', b'A'..=b'Z');

/// A transformation mapping for use with [`retain_by_recoding`] which
/// uppercases nucleotides and filters `\n\r\t :.~-`.
///
/// [`retain_by_recoding`]: RetainSequence::retain_by_recoding
const NT_HASH_FILTER: ByteMap = AA_HASH_FILTER.map_to_one(b"~", 0);

/// Returns SHA1 hex of uppercased sequence with `\n\r\t :.~-` removed.
#[inline]
pub fn nt_id(seq: &Nucleotides) -> Option<String> {
    let cleaned = filter_nucleotides(seq);
    (!cleaned.is_empty()).then(|| sha1_hex(&cleaned))
}

/// Returns MD5 hex of uppercased sequence with `\n\r\t :.-` removed (keeps
/// `~`).
#[inline]
pub fn variant_hash(seq: &AminoAcids) -> Option<String> {
    let cleaned = filter_amino_acids(seq);
    (!cleaned.is_empty()).then(|| md5_hex(&cleaned))
}

/// Filter and uppercase nucleotide sequence for hashing.
#[inline]
fn filter_nucleotides(seq: &Nucleotides) -> Vec<u8> {
    let mut bytes = seq.as_bytes().to_vec();
    bytes.retain_by_recoding(&NT_HASH_FILTER);
    bytes
}

/// Filter and uppercase amino acid sequence for hashing.
#[inline]
fn filter_amino_acids(seq: &AminoAcids) -> Vec<u8> {
    let mut bytes = seq.as_bytes().to_vec();
    bytes.retain_by_recoding(&AA_HASH_FILTER);
    bytes
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
        assert!(nt_id(&seq).is_none());
    }

    #[test]
    fn test_nt_id_only_whitespace() {
        let seq: Nucleotides = b"\n\r\t ".into();
        assert!(nt_id(&seq).is_none());
    }

    #[test]
    fn test_nt_id_only_gaps() {
        let seq: Nucleotides = b":.~-".into();
        assert!(nt_id(&seq).is_none());
    }

    #[test]
    fn test_nt_id_case_normalization() {
        let upper: Nucleotides = b"ACGT".into();
        let lower: Nucleotides = b"acgt".into();
        let mixed: Nucleotides = b"AcGt".into();

        assert_eq!(nt_id(&upper), nt_id(&lower));
        assert_eq!(nt_id(&upper), nt_id(&mixed));
    }

    #[test]
    fn test_nt_id_whitespace_removal() {
        let clean: Nucleotides = b"ACGT".into();
        let messy: Nucleotides = b"A\nC\rG\tT ".into();

        assert_eq!(nt_id(&clean), nt_id(&messy));
    }

    #[test]
    fn test_nt_id_gap_removal() {
        let clean: Nucleotides = b"ACGT".into();
        let gapped: Nucleotides = b"A:C.G~T-".into();

        assert_eq!(nt_id(&clean), nt_id(&gapped));
    }

    #[test]
    fn test_variant_hash_basic() {
        let seq: AminoAcids = b"ACDE".into();
        let hash = variant_hash(&seq).unwrap();

        assert_eq!(hash.len(), 32);
        assert!(hash.chars().all(|c| c.is_ascii_hexdigit()));
    }

    #[test]
    fn test_variant_hash_empty() {
        let seq = AminoAcids::new();
        assert!(variant_hash(&seq).is_none());
    }

    #[test]
    fn test_variant_hash_case_normalization() {
        let upper: AminoAcids = b"ACDE".into();
        let lower: AminoAcids = b"acde".into();

        assert_eq!(variant_hash(&upper), variant_hash(&lower));
    }

    #[test]
    fn test_variant_hash_whitespace_removal() {
        let clean: AminoAcids = b"ACDE".into();
        let messy: AminoAcids = b"A\nC\rD\tE ".into();

        assert_eq!(variant_hash(&clean), variant_hash(&messy));
    }

    #[test]
    fn test_variant_hash_preserves_tilde() {
        let with_tilde: AminoAcids = b"ACDE~".into();
        let without_tilde: AminoAcids = b"ACDE".into();

        assert_ne!(variant_hash(&with_tilde), variant_hash(&without_tilde));
    }

    #[test]
    fn test_variant_hash_removes_gaps() {
        let clean: AminoAcids = b"ACDE".into();
        let gapped: AminoAcids = b"A:C.D-E".into();

        assert_eq!(variant_hash(&clean), variant_hash(&gapped));
    }
}
