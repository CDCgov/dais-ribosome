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
    /// The ID of the query.
    pub id:          String,
    /// The nucleotides sequence, containing unaligned, uppercase IUPAC. `U` is
    /// preserved in addition to `T`.
    pub nucleotides: Nucleotides,
    /// The compound type of the query.
    pub ctype:       String,
}
