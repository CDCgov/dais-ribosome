use crate::conf::TSVReader;
use std::{collections::HashMap, path::Path};
use zoe::{
    data::{fasta::FastaSeq, nucleotides::ToDNA, records::RecordReader},
    prelude::*,
};

#[derive(Debug)]
pub struct BothStrands {
    pub forward: Nucleotides,
    pub reverse: Nucleotides,
}

#[derive(Hash, PartialEq, Eq, Debug)]
pub struct RefKey {
    pub reference_id:  String,
    pub compound_type: String,
}

impl From<(&str, &str)> for RefKey {
    fn from((ref_id, ctype): (&str, &str)) -> Self {
        RefKey {
            reference_id:  ref_id.to_string(),
            compound_type: ctype.to_string(),
        }
    }
}

/// Reads FASTA references from a file where each sequence name is pipe-delimited
/// in the format: `reference_id|compound_type`
pub fn get_fasta_refs(path: &Path) -> Result<HashMap<RefKey, Vec<BothStrands>>, std::io::Error> {
    let data = FastaReader::from_filename(path)?;
    let mut refs = HashMap::new();

    for r in data {
        let FastaSeq { name, sequence } = r?;

        let forward = sequence.recode_to_dna();
        let reverse = forward.to_reverse_complement();
        let strands = BothStrands { forward, reverse };

        // Parse pipe-delimited name: reference_id|compound_type
        let parts: Vec<&str> = name.split('|').collect();
        if parts.len() < 2 {
            return Err(TSVReader::new_invalid_err(
                &format!("Reference FASTA name must have format 'reference_id|compound_type', but found '{name}'",),
                path,
            ));
        }

        let key = (parts[0], parts[1]).into();
        refs.entry(key).or_insert_with(Vec::new).push(strands);
    }

    Ok(refs)
}
