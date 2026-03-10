pub use crate::annotation::error::RibosomeError;

use zoe::{
    data::{fasta::FastaSeq, nucleotides::ToDNA},
    prelude::*,
};

/// [`QueryRecord`] contains the id, compound type (ctype), and [`Nucleotides`]
/// data.
///
/// The nucleotides are guaranteed to be unaligned, uppercase IUPAC (with `T`
/// instead of `U`).
#[derive(Debug)]
pub struct QueryRecord {
    pub id:          String,
    pub nucleotides: Nucleotides,
    pub ctype:       String,
}

impl TryFrom<FastaSeq> for QueryRecord {
    type Error = RibosomeError;

    fn try_from(fasta: FastaSeq) -> Result<Self, Self::Error> {
        let FastaSeq { name, sequence } = fasta;
        // BREAKING: we previously only removed: '*: .~-'
        let nucleotides = sequence.filter_to_dna_unaligned();

        if nucleotides.is_empty() {
            return Err(RibosomeError::InvalidSequence(name.to_string()));
        }

        if name.contains('|') {
            let mut parts = name.split('|').map(|part| part.trim_ascii().to_string());

            let (Some(id), Some(ctype)) = (parts.next(), parts.next()) else {
                return Err(RibosomeError::InvalidFastaFormat);
            };

            if ctype.is_empty() {
                return Err(RibosomeError::NoCtype(name));
            }

            Ok(QueryRecord { id, nucleotides, ctype })
        } else {
            // TODO : handle unclassified queries
            Ok(QueryRecord {
                id: name.trim_ascii().to_string(),
                nucleotides,
                ctype: String::new(),
            })
        }
    }
}
