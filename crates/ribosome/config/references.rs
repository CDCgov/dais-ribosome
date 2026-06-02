//! The module data structure and parsing for it.

use crate::data::keys::RefKey;
use std::{collections::HashMap, path::Path};
use zoe::{
    data::{fasta::FastaSeq, nucleotides::ToDNA},
    prelude::{FastaReader, Nucleotides},
};

/// Loads reference sequences from a FASTA file.
///
/// Each record header must be pipe-delimited of the form:
/// `reference_id|compound_type`. Any additional pipe-delimited fields are
/// ignored.
///
/// This function also recodes the sequence to uppercase IUPAC with corrected
/// gaps, using `N` for anything that cannot be recoded.
///
/// Each `Vec` in the output [`HashMap`] is guaranteed to be non-empty.
///
/// ## Errors
///
/// All IO errors are propagated without path context. An error is also returned
/// if a sequence name doesn't match the expected format, with context including
/// the expected format and header.
pub fn load_references(path: &Path) -> std::io::Result<HashMap<RefKey, Vec<Nucleotides>>> {
    let data = FastaReader::from_path(path)?;
    let mut refs = HashMap::new();

    for r in data {
        let FastaSeq { name, sequence } = r?;

        let forward = sequence.recode_to_dna();

        let key = RefKey::parse(&name).ok_or_else(|| {
            std::io::Error::other(format!(
                "Reference FASTA header must have format '<reference_id>|<compound_type>', but found '{name}'",
            ))
        })?;

        refs.entry(key).or_insert_with(Vec::new).push(forward);
    }

    // Validity: each value in the output HashMap will be non-empty since we
    // only insert entries when there is a sequence to push
    Ok(refs)
}
