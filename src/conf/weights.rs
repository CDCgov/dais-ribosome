use crate::conf::{SpecKey, TSVReader};
use std::{
    collections::HashMap,
    fmt::Debug,
    io::{BufRead, BufReader},
    path::Path,
};
use zoe::data::records::RecordReader;

#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub struct CodonKey {
    pub position: u32,
    pub codon:    [u8; 3],
}

impl Debug for CodonKey {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{p}\t{c1}{c2}{c3}",
            p = self.position,
            c1 = self.codon[0] as char,
            c2 = self.codon[1] as char,
            c3 = self.codon[2] as char,
        )
    }
}

type CodonPositionWeights = HashMap<CodonKey, u32>;

pub fn get_codon_weight_matrix(path: &Path) -> Result<HashMap<SpecKey, CodonPositionWeights>, std::io::Error> {
    let file = TSVReader::open_nonempty_file(path)?;
    let reader = BufReader::new(file);
    let mut codon_weight_matrix = HashMap::new();
    let mut all_matrices = HashMap::new();
    let mut last_key = ("", "").into();

    for line in reader.lines().map_while(Result::ok) {
        // Start new section
        if line.contains('|') {
            let mut parts = line.split('|');
            let key: SpecKey = parts
                .next()
                .zip(parts.next())
                .ok_or_else(|| TSVReader::new_invalid_err("Issue parsing weight matrix header.", path))?
                .into();

            if !codon_weight_matrix.is_empty() {
                all_matrices.insert(last_key, codon_weight_matrix.clone());
            }
            codon_weight_matrix.clear();
            last_key = key;
        }

        let mut parts = line.split('\t');
        if let (Some(position_str), Some(codon_str), Some(count_str)) = (parts.next(), parts.next(), parts.next())
            && let (Ok(position), Ok(count), Ok(codon)) = (
                position_str.trim_ascii().parse(),
                count_str.trim_ascii().parse::<u32>(),
                codon_str.trim_ascii().as_bytes()[..3].try_into(),
            )
        {
            let mut key = CodonKey { position, codon };
            key.codon.make_ascii_uppercase();
            codon_weight_matrix.insert(key, count);
        }
    }

    if !codon_weight_matrix.is_empty() {
        all_matrices.insert(last_key, codon_weight_matrix);
    }

    Ok(all_matrices)
}
