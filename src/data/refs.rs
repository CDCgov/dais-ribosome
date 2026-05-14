//! Reference sequence loading from FASTA files.

use crate::data::keys::RefKey;
use std::{collections::HashMap, io::Error as IOError, path::Path};
use zoe::{
    data::{fasta::FastaSeq, nucleotides::ToDNA},
    prelude::*,
};

/// Map from reference key to list of reference sequences.
pub type ReferenceMap = HashMap<RefKey, Vec<Nucleotides>>;

/// Load reference sequences from a FASTA file.
///
/// Each sequence name must be pipe-delimited: `reference_id|compound_type`
///
/// # Errors
///
/// Returns an error if:
/// - The file cannot be read
/// - A sequence name doesn't match the expected format
pub fn load_references(path: &Path) -> Result<ReferenceMap, std::io::Error> {
    let data = FastaReader::from_path(path)?;
    let mut refs = HashMap::new();

    for r in data {
        let FastaSeq { name, sequence } = r?;

        let forward = sequence.recode_to_dna();

        let key = RefKey::parse(&name).ok_or_else(|| {
            IOError::other(format!(
                "Reference FASTA header must have format '<reference_id>|<compound_type>', but found '{name}'",
            ))
        })?;

        refs.entry(key).or_insert_with(Vec::new).push(forward);
    }

    Ok(refs)
}
