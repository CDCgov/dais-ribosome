//! Miscellaneous data types used by DAIS-ribosome.

pub(crate) mod coords;
pub(crate) mod exons;
pub(crate) mod keys;
pub mod ranges;
pub(crate) mod weights;

use std::{error::Error, fmt::Display};
use zoe::{
    data::{ByteMap, RetainSequence, err::GetCode},
    prelude::{Len, Nucleotides},
};

/// [`QueryRecord`] contains the id, compound type (ctype), and [`Nucleotides`]
/// data.
#[derive(Debug)]
pub struct QueryRecord {
    /// The ID of the query.
    pub(crate) id:          String,
    /// The nucleotides sequence, containing unaligned, uppercase IUPAC. `U` is
    /// preserved in addition to `T`.
    pub(crate) nucleotides: Nucleotides,
    /// The compound type of the query.
    pub(crate) ctype:       String,
}

impl QueryRecord {
    /// Forms a new [`QueryRecord`] from an `id`, a nucleotide `sequence`, and a
    /// compound type.
    ///
    /// The `sequence` is sanitized as follows: TODO
    ///
    /// ## Errors
    ///
    /// If the sequence is empty after sanitization, [`NoNucleotides`] is
    /// returned.
    pub fn new(id: String, sequence: Vec<u8>, ctype: String) -> Result<Self, NoNucleotides> {
        // BREAKING: we previously only removed: '*: .~-'
        let nucleotides = sanitize_seq(sequence);

        if nucleotides.is_empty() {
            return Err(NoNucleotides { id });
        }

        Ok(Self { id, nucleotides, ctype })
    }

    /// Returns the ID of the query.
    #[inline]
    pub fn id(&self) -> &String {
        &self.id
    }

    /// Returns the sanitized query sequence.
    #[inline]
    pub fn nucleotides(&self) -> &Nucleotides {
        &self.nucleotides
    }

    /// Returns the compound type of the query.
    #[inline]
    pub fn ctype(&self) -> &String {
        &self.ctype
    }
}

/// An error caused by the query sequence being empty after sanitization.
#[derive(Clone, Eq, PartialEq, Hash, Debug)]
pub struct NoNucleotides {
    /// The ID of the query.
    pub id: String,
}

impl Display for NoNucleotides {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "A sequence contained no unaligned DNA data. See id: {id}", id = self.id)
    }
}

impl Error for NoNucleotides {}
impl GetCode for NoNucleotides {}

/// Sanitizes an incoming sequence so that it meets the validity requirements of
/// [`QueryRecord`].
///
/// This converts to uppercase, preserves IUPAC characters, preserves `U` in
/// addition to `T`, and preserves `X`. All other bytes are removed.
#[must_use]
#[cfg(feature = "regression-testing")]
fn sanitize_seq(mut seq: Vec<u8>) -> Nucleotides {
    const SANITIZE: ByteMap = ByteMap::all(0)
        .preserve_range(b'A'..=b'Z')
        .preserve_range(b'a'..=b'z')
        .map(b"acgturyswkmbdhvn", b"ACGTURYSWKMBDHVN")
        .map(b"ux", b"UX");

    seq.retain_by_recoding(&SANITIZE);
    Nucleotides::from(seq)
}

/// Sanitizes an incoming sequence so that it meets the validity requirements of
/// [`QueryRecord`].
///
/// This converts to uppercase, preserves IUPAC characters, preserves `U` in
/// addition to `T`, and maps `X` to `N`. All other bytes are removed.
#[must_use]
#[cfg(not(feature = "regression-testing"))]
fn sanitize_seq(mut seq: Vec<u8>) -> Nucleotides {
    const SANITIZE: ByteMap = ByteMap::all(0)
        .preserve_range(b'A'..=b'Z')
        .preserve_range(b'a'..=b'z')
        .map(b"acgturyswkmbdhvn", b"ACGTURYSWKMBDHVN")
        .map(b"uxX", b"UNN");

    seq.retain_by_recoding(&SANITIZE);
    Nucleotides::from(seq)
}
