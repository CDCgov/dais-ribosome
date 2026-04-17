use crate::{
    data::{
        exons::{ExonCoords, Exons},
        keys::RefKey,
    },
    ranges::RangeExt,
};
use std::{
    cmp::Ordering,
    collections::HashMap,
    fs::File,
    io::{BufRead, BufReader, ErrorKind, Lines},
    iter::Enumerate,
    ops::Range,
    path::Path,
    str::FromStr,
};
use zoe::{
    data::{
        SanitizeBase,
        err::{ResultWithErrorContext, WithErrorContext},
    },
    prelude::RefineDNAStrat,
    unwrap_or_return_some_err,
};

/// A hash map from [`RefKey`] values to a vector of the proteins (e.g., `HA`,
/// `HA-signal`) and their [`Exons`].
pub(crate) type CdsSpecMap = HashMap<RefKey, Vec<(String, Exons)>>;

/// Loads the coding sequence specifications from a TSV file.
///
/// The expected columns are `reference_id`, `protein`, `ctype`, `coords`, and
/// then optionally `required_start`.
///
/// The coordinates are in `start..end` format (1-based inclusive-end), with
/// multiple exons separated by semicolons.
///
/// ## Errors
///
/// Returns an error (without path context) if:
///
/// - The file cannot be read
/// - A parsing error occurs on one of the lines (e.g., missing required field,
///   invalid range, invalid residues in required beginning field, etc.)
pub(crate) fn load_cds_spec(path: &Path) -> std::io::Result<CdsSpecMap> {
    let reader = TsvReader::from_path(path)?;

    let mut cds_specs: HashMap<RefKey, Vec<(String, Exons)>> = HashMap::new();

    for row in reader {
        let TsvRow {
            reference_id,
            protein,
            ctype,
            coords,
            mut required_start,
        } = row?;

        // Standardize the required start
        if let Some(codon) = &mut required_start {
            for residue in codon {
                let Some(new_residue) = residue.refine_base(RefineDNAStrat::AcgtNoGapsUc) else {
                    return Err(std::io::Error::other(format!(
                        "An invalid residue ({residue}) was found in the required beginning field. Expected any of ACGTUacgtu",
                        residue = *residue as char
                    )));
                };
                *residue = new_residue;
            }
        }

        let cds_len = coords
            .last()
            .expect("The coords field of TsvRow should be non-empty")
            .cds_range
            .end;

        if !cds_len.is_multiple_of(3) {
            return Err(std::io::Error::other(
                "The length of the coding sequence (sum of all exon lengths) was not a multiple of 3.",
            )
            .with_context(format!("Failed to parse {ctype} product for reference ID {reference_id}"))
            .into());
        }

        // Validity: coords is non-empty, CDS length is multiple of 3, etc.
        let exons = Exons { required_start, coords };

        let key = RefKey::new(reference_id, ctype);
        cds_specs.entry(key).or_default().push((protein, exons));
    }

    Ok(cds_specs)
}

// TODO: Switch this to returning ModuleLoadError directly?

/// The parsed data from a single row of the `cds-spec.tsv` file for the module.
#[derive(Clone, Debug)]
struct TsvRow {
    /// The reference ID in column 1 (e.g., `ANHUI01`, `PHUKET3073`).
    reference_id: String,

    /// The compound type in column 2 (e.g., `A_HA_H7`, `B_HA`).
    ctype: String,

    /// The protein product name in column 3 (e.g., `HA`, `HA-signal`).
    protein: String,

    /// A list of the exon coordinates in column 4 (e.g., `55..1683`,
    /// `1..33;689..1024`).
    ///
    /// The exons are ordered by `cds_range`, which form a partition of
    /// `0..cds_len` where `cds_len` is the total length of the coding sequence.
    /// The `ref_range` fields are in order, although up to 2 nucleotides
    /// overlap is allowed between ranges. Note that any repeated indices are
    /// represented twice with distinct coordinates in the coding sequence.
    ///
    /// This vector is non-empty.
    coords: Vec<ExonCoords>,

    /// An optional required codon at the start in column 5 (e.g., `ATG`).
    required_start: Option<[u8; 3]>,
}

impl FromStr for TsvRow {
    type Err = std::io::Error;

    /// Parses a string `line` to a [`TsvRow`]. The string must be trimmed and
    /// non-empty, and should not be the header row or a commented line.
    fn from_str(line: &str) -> std::io::Result<Self> {
        // End the iterator on empty fields, so that missing field errors appear
        let mut parts = line.split('\t').map(str::trim_ascii).take_while(|s| !s.is_empty());

        // Get fields as string slices
        let Some(reference_id) = parts.next() else {
            return Err(std::io::Error::other("Missing Reference ID field (first field)"));
        };

        let Some(ctype) = parts.next() else {
            return Err(std::io::Error::other("Missing Compound Type field (second field)"));
        };

        let Some(protein) = parts.next() else {
            return Err(std::io::Error::other("Missing Protein field (third field)"));
        };

        let Some(coords) = parts.next() else {
            return Err(std::io::Error::other("Missing Coords field (fourth field)"));
        };
        let required_start = parts.next();

        let reference_id = reference_id.to_string();
        let protein = protein.to_string();
        let ctype = ctype.to_string();

        // Parse coordinate ranges (e.g., 1..54 or 1..36;692..1027)
        let coords = parse_coordinate_ranges(coords)?;

        // Parse optional required start codon
        let required_start = match required_start {
            Some(required_start) => {
                let mut codon: [u8; 3] = required_start.as_bytes().try_into().map_err(|_| {
                    std::io::Error::other(format!(
                        "If provided, Required Beginning must contain exactly 3 amino acids (not {})",
                        required_start.len()
                    ))
                })?;
                codon.make_ascii_uppercase();
                Some(codon)
            }
            None => None,
        };

        Ok(TsvRow {
            reference_id,
            protein,
            ctype,
            coords,
            required_start,
        })
    }
}

/// A reader over the `cds-spec.tsv` file of a module, automatically parsing
/// each line into a [`TsvRow`].
///
/// ## Errors
///
/// Any IO errors are propagated without context. Context with the line number
/// is added for parsing errors.
struct TsvReader {
    lines: Enumerate<Lines<BufReader<File>>>,
}

impl TsvReader {
    fn from_path(path: impl AsRef<Path>) -> std::io::Result<Self> {
        let path = path.as_ref();
        let file = File::open(path).with_path_context("Failed to open path", path)?;
        let mut reader = std::io::BufReader::new(file);

        if reader.fill_buf()?.is_empty() {
            return Err(std::io::Error::new(
                ErrorKind::InvalidData,
                format!("No data was found at path {path}", path = path.display()),
            ));
        }

        Ok(Self {
            lines: reader.lines().enumerate(),
        })
    }
}

impl Iterator for TsvReader {
    type Item = std::io::Result<TsvRow>;

    fn next(&mut self) -> Option<Self::Item> {
        let (line_idx, line) = self.lines.next()?;
        let line = unwrap_or_return_some_err!(line);
        let line = line.trim();

        // Skip empty lines and headers
        if line.is_empty() || line.starts_with("Reference ID") || line.starts_with('#') {
            return self.next();
        }

        Some(
            line.parse()
                .with_context(format!("Failed to parse line {line_num}: {line}", line_num = line_idx + 1))
                .map_err(Into::into),
        )
    }
}

/// Parses a non-empty semicolon-delimited list of 1-based end-inclusive
/// reference ranges, converting them to 0-based [`ExonCoords`] ranges.
///
/// The length of the coding sequence is the sum of the lengths of all the
/// reference ranges. The `cds_range` fields of the resulting [`ExonCoords`]
/// partition this length, starting from 0.
///
/// The output vector will be non-empty.
///
/// ## Errors
///
/// - Each range must successfully parse
/// - The ranges must be in order, with at most 2 nt of overlap
fn parse_coordinate_ranges(coords: &str) -> std::io::Result<Vec<ExonCoords>> {
    /// The maximum amount of overlap allowed between exons.
    ///
    /// SARS-CoV-2 requires -1 exon-to-exon frameshift with other viruses
    /// reported up to -2.
    const MAX_DUPLICATED_OVERLAP_NT: usize = 2;

    let mut exon_ranges: Vec<ExonCoords> = Vec::new();
    let mut cds_start = 0;

    // Parses and pushes a range to exon_ranges. Iterators do not work since we
    // need to access the previous ExonCoords to validate order and overlap.
    for coord_range in coords.split(';') {
        let ref_range = parse_coordinate_range(coord_range)?;

        if let Some(last) = exon_ranges.last() {
            match ref_range.relaxed_cmp(&last.ref_range) {
                Some(Ordering::Greater) => {}
                Some(Ordering::Less) => {
                    return Err(std::io::Error::other(format!(
                        "Exons out of order! Found {} then {}",
                        last.ref_range.display_inclusive(),
                        ref_range.display_inclusive(),
                    )));
                }
                Some(Ordering::Equal) => {
                    return Err(std::io::Error::other(format!(
                        "Found the same exon twice: {}",
                        ref_range.display_inclusive()
                    )));
                }
                None => {
                    return Err(std::io::Error::other(format!(
                        "One exon cannot completely contain another! Found {} then {}",
                        last.ref_range.display_inclusive(),
                        ref_range.display_inclusive(),
                    )));
                }
            }

            // Exclusive index - inclusive index is valid length
            let overlap_nt = last.ref_range.end.saturating_sub(ref_range.start);
            if overlap_nt > MAX_DUPLICATED_OVERLAP_NT {
                return Err(std::io::Error::other(format!(
                    "Exon overlap exceeds {MAX_DUPLICATED_OVERLAP_NT} nt! Found {} then {}",
                    last.ref_range.display_inclusive(),
                    ref_range.display_inclusive(),
                )));
            }
        }

        let cds_end = cds_start + ref_range.len();

        // Validity: ref_range is non-empty per guarantees from
        // parse_coordinate_range, and they are the same length per above
        // definition
        exon_ranges.push(ExonCoords {
            ref_range,
            cds_range: cds_start..cds_end,
        });
        cds_start = cds_end;
    }

    // Validity: exon_ranges will be non-empty since split is always non-empty.
    Ok(exon_ranges)
}

/// Parses a string containing a 1-based end-inclusive range, converting it to a
/// 0-based [`Range`].
///
/// The returned range is guaranteed to be non-empty.
fn parse_coordinate_range(coord_range: &str) -> std::io::Result<Range<usize>> {
    let coord_range = coord_range.trim();
    let range_parts: Vec<&str> = coord_range.split("..").collect();

    let &[start_part, end_part] = range_parts.as_slice() else {
        return Err(std::io::Error::other(format!("Invalid coordinate range '{coord_range}'")));
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
            "End coordinate must be >= start ({start}..{end})",
        )));
    }

    // Convert to 0-based half-open range (inclusive start, exclusive end)
    let Some(start) = start.checked_sub(1) else {
        return Err(std::io::Error::other("Start coordinate must be at least 1"));
    };

    Ok(start..end)
}
