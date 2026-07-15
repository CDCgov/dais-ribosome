//! The module data structure and parsing for it.

use crate::data::keys::RefKey;
use std::{collections::HashMap, path::Path};
use zoe::{
    data::{
        err::{ResultWithErrorContext, WithErrorContext},
        fasta::FastaSeq,
        nucleotides::CheckNucleotides,
    },
    prelude::{FastaReader, IsValidDNA, Nucleotides},
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
/// - All IO errors are propagated without path context.
/// - The reference ID and compound type cannot contain tabs.
/// - The sequence name must match the expected format.
/// - The sequence must contain IUPAC bases, and no gaps can be present.
pub fn load_references(path: &Path) -> std::io::Result<HashMap<RefKey, Vec<Nucleotides>>> {
    let data = FastaReader::from_path(path)?;
    let mut refs = HashMap::new();

    for r in data {
        let FastaSeq { name, sequence } = r?;

        let sequence = Nucleotides::from(sequence);

        if !sequence.is_valid_dna(IsValidDNA::IupacNoGaps) {
            let e =
                if let Some(p) = zoe::search::position_by_byte2::<32>(sequence.as_bytes(), b'-', b'.') {
                    std::io::Error::other(format!(
                        "A gap character '{}' was found. Consider N-filling gaps instead.",
                        sequence[p] as char
                    ))
                } else {
                    std::io::Error::other(
                        "Non-IUPAC character(s) were found. Consider removing them (if it is erroneously present) or masking to N.",
                    )
                }
            .with_context(format!("Failed to parse reference sequence with name: {name}"));

            return Err(e.into());
        }

        let key = RefKey::parse(&name).with_context("Failed to parse FASTA header")?;

        refs.entry(key).or_insert_with(Vec::new).push(sequence);
    }

    // Validity: each value in the output HashMap will be non-empty since we
    // only insert entries when there is a sequence to push
    Ok(refs)
}
