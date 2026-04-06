//! Miscellaneous data types used by DAIS-ribosome.

pub(crate) mod coords;
pub(crate) mod exons;
pub(crate) mod keys;
pub mod ranges;
pub(crate) mod weights;

use zoe::prelude::Nucleotides;

/// [`QueryRecord`] contains the id, compound type (ctype), and [`Nucleotides`]
/// data.
#[derive(Debug)]
pub struct QueryRecord {
    pub id:          String,
    /// The nucleotides sequence, containing unaligned, uppercase IUPAC (with `T`
    /// instead of `U`).
    pub nucleotides: Nucleotides,
    pub ctype:       String,
}
