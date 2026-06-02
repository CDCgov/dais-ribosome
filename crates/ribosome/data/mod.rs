//! Miscellaneous data types used by DAIS-ribosome.

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
    id:          String,
    /// The nucleotides sequence, containing unaligned, uppercase IUPAC. `U` is
    /// preserved in addition to `T`.
    nucleotides: Nucleotides,
    /// The compound type of the query.
    ctype:       String,
}

impl QueryRecord {
    /// Forms a new [`QueryRecord`] from an `id`, a nucleotide `sequence`, and a
    /// compound type.
    ///
    /// The sequence is sanitized as follows:
    ///
    /// - DNA IUPAC is preserved/uppercased (including preserving `U`)
    /// - Any other alphabetic bytes are mapped to uppercase `N`
    /// - Any non-alphabetic bytes are filtered
    ///
    /// ## Errors
    ///
    /// If the sequence is empty after sanitization, [`NoNucleotides`] is
    /// returned.
    pub fn new(id: String, sequence: Vec<u8>, ctype: String) -> Result<Self, NoNucleotides> {
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

    /// Returns the owned compound type of the query.
    pub(crate) fn into_ctype(self) -> String {
        self.ctype
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
/// This preserves/uppercases DNA IUPAC (including preserving `U`), recodes any
/// other alphabetic bytes to uppercase `N`, and filters any non-alphabetic
/// bytes.
#[must_use]
fn sanitize_seq(mut seq: Vec<u8>) -> Nucleotides {
    const SANITIZE: ByteMap = ByteMap::all(0)
        .map_range_to_one(b'A'..=b'Z', b'N')
        .map_range_to_one(b'a'..=b'z', b'N')
        .preserve(b"ACGTURYSWKMBDHVN")
        .map(b"acgturyswkmbdhvn", b"ACGTURYSWKMBDHVN");

    seq.retain_by_recoding(&SANITIZE);
    Nucleotides::from(seq)
}
