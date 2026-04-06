pub(crate) mod coords;
pub(crate) mod ctype;
pub(crate) mod exons;
pub(crate) mod keys;
mod outputs;
pub(crate) mod products;
pub(crate) mod ranges;
pub(crate) mod refs;
pub(crate) mod spec;
pub(crate) mod weights;

pub use outputs::*;

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
