use crate::conf::TSVReader;
use std::{collections::HashMap, io::BufRead, path::Path};
use zoe::data::records::RecordReader;

#[derive(Debug)]
pub struct Exons {
    pub ctype:          String,
    pub required_start: Option<[u8; 3]>,
    pub ranges:         Vec<std::ops::Range<usize>>,
}

#[derive(Hash, PartialEq, Eq, Debug)]
pub struct SpecKey {
    pub reference_id:    String,
    pub protein_product: String,
}

impl From<(&str, &str)> for SpecKey {
    fn from((ref_id, prot): (&str, &str)) -> Self {
        SpecKey {
            reference_id:    ref_id.to_string(),
            protein_product: prot.to_string(),
        }
    }
}

pub fn get_cds_spec(path: &Path) -> Result<HashMap<SpecKey, Exons>, std::io::Error> {
    let file = TSVReader::open_nonempty_file(path)?;
    let reader = std::io::BufReader::new(file);
    let mut result = HashMap::new();

    for line in reader.lines() {
        let line = line.map_err(|err| TSVReader::new_wrapped("line read failed", path, err))?;
        let line = line.trim();

        if line.is_empty() || line.starts_with("Reference ID") || line.starts_with('#') {
            continue;
        }

        let mut parts = line.split('\t');
        let (Some(reference_id), Some(ctype), Some(protein), Some(coords)) =
            (parts.next(), parts.next(), parts.next(), parts.next())
        else {
            return Err(TSVReader::new_invalid_err("Missing required field(s)", path));
        };

        let reference_id = reference_id.trim_ascii();
        let protein = protein.trim_ascii();
        let ctype = ctype.trim_ascii().to_string();
        let coords = coords.trim_ascii().trim_matches('"');

        let required_start = parts
            .next()
            .map(|s| s.trim_ascii().as_bytes())
            .and_then(|s| s.split_first_chunk::<3>())
            .map(|ch| {
                let mut c = *ch.0;
                c.make_ascii_uppercase();
                c
            });

        // Parse coordinate ranges (e.g., "1..54" or "1..36;692..1027")
        let mut ranges = Vec::new();
        for coord_range in coords.split(';') {
            let coord_range = coord_range.trim();
            let range_parts: Vec<&str> = coord_range.split("..").collect();
            if range_parts.len() != 2 {
                return Err(TSVReader::new_invalid_err(
                    &format!("Invalid coordinate range '{coord_range:?}'"),
                    path,
                ));
            }
            let start: usize = range_parts[0].parse().map_err(|_| {
                TSVReader::new_invalid_err(&format!("Invalid start coordinate '{s}'", s = range_parts[0]), path)
            })?;

            let end: usize = range_parts[1].parse().map_err(|_| {
                TSVReader::new_invalid_err(&format!("Invalid end coordinate '{e}'", e = range_parts[1]), path)
            })?;

            if start < 1 {
                return Err(TSVReader::new_invalid_err("Start coordinate must be at least 1 (0)", path));
            }
            if end < start {
                return Err(TSVReader::new_invalid_err(
                    &format!("End coordinate must be >= start ({start}..{end})"),
                    path,
                ));
            }

            // Convert from 1-based inclusive to 0-based half-open range
            ranges.push((start - 1)..end);
        }

        let key = (reference_id, protein).into();
        result.insert(
            key,
            Exons {
                ctype,
                required_start,
                ranges,
            },
        );
    }

    Ok(result)
}
