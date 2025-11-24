//! CDS specification loading.

use crate::data::{
    TSVReader,
    exons::{CtypeExons, ExonCoords},
    keys::SpecKey,
};
use std::{
    collections::HashMap,
    io::{BufRead, Error as IOError},
    path::Path,
};
use zoe::data::err::ResultWithErrorContext;

/// Map from spec key to exon specification (with ctype for grouping).
pub type CdsSpecMap = HashMap<SpecKey, CtypeExons>;

/// Load CDS specifications from a TSV file.
///
/// Expected columns: reference_id, protein, ctype, coords, [required_start]
///
/// Coordinates are in "start..end" format (1-based inclusive), with multiple
/// exons separated by semicolons.
///
/// # Errors
///
/// Returns an error if:
/// - The file cannot be read
/// - A line is missing required fields
/// - Coordinate ranges are invalid
pub fn load_cds_spec(path: &Path) -> Result<CdsSpecMap, std::io::Error> {
    let file = TSVReader::open_nonempty_file(path)?;
    let reader = std::io::BufReader::new(file);
    let mut result = HashMap::new();

    for line in reader.lines() {
        let line = line.with_file_context("line read failed", path)?;
        let line = line.trim();

        // Skip empty lines and headers
        if line.is_empty() || line.starts_with("Reference ID") || line.starts_with('#') {
            continue;
        }

        let mut parts = line.split('\t');
        let (Some(reference_id), Some(ctype), Some(protein), Some(coords)) =
            (parts.next(), parts.next(), parts.next(), parts.next())
        else {
            return Err(IOError::other(format!("Missing required field(s): {}", path.display())));
        };

        let reference_id = reference_id.trim_ascii();
        let protein = protein.trim_ascii();
        let ctype = ctype.trim_ascii().to_string();
        let coords = coords.trim_ascii().trim_matches('"');

        // Parse optional required start codon
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
        let coords = parse_coordinate_ranges(coords, path)?;
        let key = (reference_id, protein).into();
        result.insert(
            key,
            CtypeExons {
                ctype,
                required_start,
                coords,
            },
        );
    }

    Ok(result)
}

fn parse_coordinate_ranges(coords: &str, path: &Path) -> Result<Vec<ExonCoords>, std::io::Error> {
    let mut exon_ranges: Vec<ExonCoords> = vec![];
    let mut ref_to_cds_offset = 0;

    for coord_range in coords.split(';') {
        let coord_range = coord_range.trim();
        let range_parts: Vec<&str> = coord_range.split("..").collect();

        if range_parts.len() != 2 {
            return Err(IOError::other(format!(
                "Invalid coordinate range '{coord_range}': {}",
                path.display()
            )));
        }

        let mut start: usize = range_parts[0]
            .parse()
            .with_context(format!("Invalid start coordinate '{}'", range_parts[0]))?;

        let end: usize = range_parts[1]
            .parse()
            .with_context(format!("Invalid end coordinate '{}'", range_parts[1]))?;

        if start < 1 {
            return Err(IOError::other(format!(
                "Start coordinate must be at least 1: {}",
                path.display()
            )));
        }
        if end < start {
            return Err(IOError::other(format!(
                "End coordinate must be >= start ({start}..{end}): {}",
                path.display(),
            )));
        }

        // 1-based inclusive to 0-based inclusive
        start -= 1;
        ref_to_cds_offset += if let Some(last) = exon_ranges.iter().last() {
            // 0-based inclusive - exclusive
            start - last.ref_range.end
        } else {
            start
        };

        exon_ranges.push(ExonCoords {
            ref_range: start..end,
            ref_to_cds_offset,
        });
    }

    if exon_ranges.is_empty() {
        return Err(IOError::other(format!("No exon ranges found were found: {}", path.display())));
    }

    Ok(exon_ranges)
}
