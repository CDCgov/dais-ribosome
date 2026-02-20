//! CDS specification loading.

use crate::data::{
    TSVReader,
    exons::{CtypeExons, ExonCoords},
    keys::SpecKey,
};
use std::{collections::HashMap, io::BufRead, ops::Range, path::Path};
use zoe::data::err::ResultWithErrorContext;

/// Map from spec key to exon specification (with ctype for grouping).
pub type CdsSpecMap = HashMap<SpecKey, CtypeExons>;

// TODO: Switch this to returning ModuleLoadError directly?

/// Load CDS specifications from a TSV file.
///
/// Expected columns: reference_id, protein, ctype, coords, [required_start]
///
/// Coordinates are in "start..end" format (1-based inclusive), with multiple
/// exons separated by semicolons.
///
/// ## Errors
///
/// Returns an error if:
///
/// - The file cannot be read
/// - A line is missing required fields
/// - Coordinate ranges are invalid
pub fn load_cds_spec(path: &Path) -> Result<CdsSpecMap, std::io::Error> {
    // TODO: We aren't actually making any sort of TSVReader... this abstraction
    // holds no utility currently.
    let file = TSVReader::open_nonempty_file(path)?;
    let reader = std::io::BufReader::new(file);
    let mut result = HashMap::new();

    for line in reader.lines() {
        let line = line.with_path_context("line read failed", path)?;
        let line = line.trim();

        // Skip empty lines and headers
        if line.is_empty() || line.starts_with("Reference ID") || line.starts_with('#') {
            continue;
        }

        let mut parts = line.split('\t').map(str::trim_ascii);
        let (Some(reference_id), Some(ctype), Some(protein), Some(coords)) =
            (parts.next(), parts.next(), parts.next(), parts.next())
        else {
            return Err(std::io::Error::other(format!(
                "Missing required field(s): {}",
                path.display()
            )));
        };

        let ctype = ctype.to_string();

        // TODO: Why are we trimming quotation marks?
        let coords = coords.trim_matches('"');
        // Parse coordinate ranges (e.g., "1..54" or "1..36;692..1027")
        let coords = parse_coordinate_ranges(coords, path)?;

        // TODO: Error on invalid contents in Required Beginning, don't just
        // silently ignore. And if it is longer than a single codon, perhaps
        // also error.

        // TODO: Can we use AminoAcidsView to make any of this more clear?

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

        let key = SpecKey::new(reference_id, protein);
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

/// Parses a semicolon-delimited list of 1-based ranges, converting them to
/// 0-based [`ExonCoords`] ranges. (TODO)
fn parse_coordinate_ranges(coords: &str, path: &Path) -> Result<Vec<ExonCoords>, std::io::Error> {
    let mut exon_ranges: Vec<ExonCoords> = Vec::new();
    let mut ref_to_cds_offset = 0;

    // Parses and pushes a range to exon_ranges. Iterators do not work since we
    // need to access the previous ExonCoords to compute offset.
    for coord_range in coords.split(';') {
        let ref_range = parse_coordinate_range(coord_range, path)?;

        // Compute offset from previous exon (i.e., length of intron)
        ref_to_cds_offset += if let Some(last) = exon_ranges.last() {
            // 0-based inclusive minus exclusive will yield length of intron
            ref_range.start - last.ref_range.end
        } else {
            ref_range.start
        };

        exon_ranges.push(ExonCoords {
            ref_range,
            ref_to_cds_offset,
        });
    }

    // Validity: exon_ranges will be non-empty since split is always non-empty.
    // An empty field manifests as an error in parse_coordinate_range
    Ok(exon_ranges)
}

/// Parses a string containing a 1-based inclusive range, converting it to a
/// 0-based [`Range`].
fn parse_coordinate_range(coord_range: &str, path: &Path) -> std::io::Result<Range<usize>> {
    let coord_range = coord_range.trim();
    let range_parts: Vec<&str> = coord_range.split("..").collect();

    let &[start_part, end_part] = range_parts.as_slice() else {
        // TODO: Have path context added where this function is called
        return Err(std::io::Error::other(format!(
            "Invalid coordinate range '{coord_range}': {}",
            path.display()
        )));
    };

    // Parse 1-based inclusive range
    let start: usize = start_part
        .parse()
        .with_context(format!("Invalid start coordinate '{}'", range_parts[0]))?;
    let end: usize = end_part
        .parse()
        .with_context(format!("Invalid end coordinate '{}'", range_parts[1]))?;

    // Since we are using inclusive range, this also requires the range is
    // non-empty
    if end < start {
        return Err(std::io::Error::other(format!(
            "End coordinate must be >= start ({start}..{end}): {}",
            path.display(),
        )));
    }

    // Convert to 0-based half-open range (inclusive start, exclusive end)
    let Some(start) = start.checked_sub(1) else {
        return Err(std::io::Error::other(format!(
            "Start coordinate must be at least 1: {}",
            path.display()
        )));
    };

    Ok(start..end)
}
