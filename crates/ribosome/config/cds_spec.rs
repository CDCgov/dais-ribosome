use crate::{
    data::{
        exons::{ExonCoords, ExonOverlapCoords, Exons, NoncodingCoords},
        keys::RefKey,
    },
    ranges::{InclusiveDisplay, RangeExt, parse_coords_inclusive},
};
use std::{
    cmp::Ordering,
    collections::HashMap,
    fs::File,
    io::{BufRead, BufReader, Lines, Read},
    iter::Enumerate,
    ops::Range,
    path::Path,
    str::FromStr,
};
use zoe::{
    data::{
        SanitizeBase,
        err::{ErrorWithContext, ResultWithErrorContext, WithErrorContext, WithSubitem},
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
/// All IO errors are propagated without path context. Each row of the file must
/// successfully parse using [`TsvRow::from_str`], and there must be at least
/// one row. The length of the coding sequence (sum of all exon lengths) must be
/// a multiple of 3. The header row (if present) must contain the expected
/// columns in the correct order, although `"required_start"` is optional in the
/// header.
pub(crate) fn load_cds_spec(path: &Path) -> std::io::Result<CdsSpecMap> {
    let reader = TsvReader::from_path(path)?;

    let mut cds_specs: CdsSpecMap = HashMap::new();

    for row in reader {
        let TsvRow {
            reference_id,
            product_name,
            ctype,
            coords,
            required_start,
        } = row?;

        // Validity: coords is non-empty per TsvRow guarantees
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

        let mut overlapped_regions = Vec::new();
        let mut noncoding_regions = Vec::new();

        for [exon1, exon2] in coords.array_windows() {
            // Overlapping and containing a noncoding region between them are
            // mutually exclusive, so use "else if"
            if let Some(overlapped) = ExonOverlapCoords::new(exon1, exon2) {
                overlapped_regions.push(overlapped);
            } else if let Some(noncoding) = NoncodingCoords::new(exon1, exon2) {
                noncoding_regions.push(noncoding);
            }
        }

        // Validity: coords is non-empty, CDS length is multiple of 3, etc.
        let exons = Exons {
            // Validity: this contains solely ACGT by TsvRow guarantees
            required_start,
            // Validity: this is non-empty by TsvRow guarantees
            coords,
            overlapped_regions,
            noncoding_regions,
        };

        let key = RefKey::new(reference_id, ctype);
        cds_specs.entry(key).or_default().push((product_name, exons));
    }

    if cds_specs.is_empty() {
        return Err(std::io::Error::other("No specifications were found in the file"));
    }

    Ok(cds_specs)
}

/// The parsed data from a single row of the `cds-spec.tsv` file for the module.
#[derive(Clone, Debug)]
struct TsvRow {
    /// The reference ID in column 1 (e.g., `ANHUI01`, `PHUKET3073`).
    reference_id: String,

    /// The compound type in column 2 (e.g., `A_HA_H7`, `B_HA`).
    ctype: String,

    /// The protein product name in column 3 (e.g., `HA`, `HA-signal`).
    product_name: String,

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
    ///
    /// The nucleotides will be in `ACGT`.
    required_start: Option<[u8; 3]>,
}

impl FromStr for TsvRow {
    type Err = std::io::Error;

    /// Parses a `line` to a [`TsvRow`]. The string must be trimmed and
    /// non-empty, and not be the header row or a commented line.
    ///
    /// ## Errors
    ///
    /// If any of the first four columns are missing or empty, then an error is
    /// returned. The exon coordinates must successfully parse using
    /// [`parse_exon_coords`]. If provided, the required first codon must
    /// contain exactly 3 nucleotides in `ACGTUacgtu`.
    fn from_str(line: &str) -> std::io::Result<Self> {
        // End the iterator on empty fields, so that missing field errors appear
        let mut parts = line.split('\t').map(str::trim_ascii).take_while(|s| !s.is_empty());

        // Get fields as string slices
        let Some(reference_id) = parts.next() else {
            return Err(std::io::Error::other("Missing reference_id field (first field)"));
        };

        let Some(ctype) = parts.next() else {
            return Err(std::io::Error::other("Missing ctype field (second field)"));
        };

        let Some(product_name) = parts.next() else {
            return Err(std::io::Error::other("Missing product_name field (third field)"));
        };

        let Some(coords) = parts.next() else {
            return Err(std::io::Error::other("Missing coords field (fourth field)"));
        };
        let required_start = parts.next();

        let reference_id = reference_id.to_string();
        let product_name = product_name.to_string();
        let ctype = ctype.to_string();

        // Parse coordinate ranges (e.g., 1..54 or 1..36;692..1027)
        let coords = parse_exon_coords(coords)?;

        // Parse optional required start codon
        let required_start = match required_start {
            Some(required_start) => {
                // Validate length of codon
                let mut codon: [u8; 3] = required_start.as_bytes().try_into().map_err(|_| {
                    std::io::Error::other(format!(
                        "If provided, required_beginning must contain exactly 3 nucleotides (not {})",
                        required_start.len()
                    ))
                })?;

                // Attempt to recode ACGTUacgtu to ACGT
                for residue in &mut codon {
                    let Some(recoded) = residue.refine_base(RefineDNAStrat::AcgtNoGapsUc) else {
                        return Err(std::io::Error::other(format!(
                            "An invalid residue ({residue}) was found in the required_beginning field. Expected any of ACGTUacgtu",
                            residue = *residue as char
                        )));
                    };

                    *residue = recoded;
                }

                Some(codon)
            }
            None => None,
        };

        Ok(TsvRow {
            reference_id,
            product_name,
            ctype,
            // Validity: coords is non-empty since parse_exon_coords ensures
            // this
            coords,
            // Validity: we recoded to ACGT above
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
    /// Opens a TSV reader for the CDS specifications file from a given path.
    ///
    /// No validation is performed to ensure the file is non-empty. Since the
    /// file is slurped, this logic is handled after parsing.
    ///
    /// ## Errors
    ///
    /// IO errors opening the file are propagated without context. The header
    /// row must be well-formed if present.
    fn from_path(path: &Path) -> std::io::Result<Self> {
        let file = File::open(path)?;
        let mut reader = std::io::BufReader::new(file);

        Self::skip_header(&mut reader)?;

        Ok(Self {
            lines: reader.lines().enumerate(),
        })
    }

    /// A helper function to pre-emptively skip past the header row, if it is
    /// present.
    ///
    /// ## Errors
    ///
    /// IO errors are propagated. If `"reference_id"`, `"ctype"`,
    /// `"product_name"`, `"coords"`, or `"required_beginning"` are present in
    /// the first line, then the first line must equal the expected header (with
    /// `"required_beginning"` optionally present or excluded).
    fn skip_header(reader: &mut BufReader<File>) -> std::io::Result<()> {
        const REQUIRED: &str = "reference_id\tctype\tproduct_name\tcoords";
        const OPTIONAL: &str = "\trequired_beginning";
        const HEADER_FIELDS: [&str; 5] = ["reference_id", "ctype", "product_name", "coords", "required_beginning"];

        // Handle case where the header appears to be fully present
        if reader.starts_with(REQUIRED.as_bytes())? {
            let mut lines = reader.lines();

            // Advance past the header using lines.next
            let Some(header_line) = lines.next().transpose()? else {
                // Should be unreachable, since the underlying buffer contains
                // HEADER
                return Ok(());
            };

            if !header_line.is_ascii() {
                return Err(std::io::Error::other("The header row contained a non-ASCII character"));
            }

            if header_line.trim_ascii().contains(' ') {
                return Err(std::io::Error::other(
                    "The header row contained a space, but only tabs were expected",
                ));
            }

            let optional = header_line.strip_prefix(REQUIRED).unwrap_or_default().trim_ascii();

            if !optional.is_empty() && optional != OPTIONAL.trim_ascii() {
                return Err(ErrorWithContext::new("Incorrect header found")
                    .with_subitem(format!("Expected header: {REQUIRED}{OPTIONAL}"))
                    .with_subitem(format!("Found header:    {header_line}"))
                    .into());
            }

            return Ok(());
        }

        // Handle case where the header seems to be present, but perhaps is
        // malformed or out of order
        let buffer = reader.fill_buf()?;

        let Some(header_line) = buffer.lines().next().transpose()? else {
            // No header present, since there are no lines
            return Ok(());
        };

        if header_line
            .split('\t')
            .map(str::trim_ascii)
            .any(|field| HEADER_FIELDS.contains(&field))
        {
            return Err(ErrorWithContext::new("Incorrect header found")
                .with_subitem(format!("Expected header: {REQUIRED}{OPTIONAL}"))
                .with_subitem(format!("Found header:    {header_line}"))
                .into());
        }

        // No header present
        Ok(())
    }
}

impl Iterator for TsvReader {
    type Item = std::io::Result<TsvRow>;

    /// Parses the next row of the CDS specifications file.
    ///
    /// Empty lines and commented lines (those beginning with `#`) are skipped.
    ///
    /// ## Errors
    ///
    /// Any errors parsing the row with [`TsvRow::from_str`] are propagated with
    /// context (including the line and line number).
    fn next(&mut self) -> Option<Self::Item> {
        for (line_idx, line) in self.lines.by_ref() {
            let line = unwrap_or_return_some_err!(line);
            let line = line.trim();

            // Ensure line is non-empty and is not a comment line
            if !line.is_empty() && !line.starts_with('#') {
                return Some(
                    line.parse()
                        .with_context(format!("Failed to parse line {line_num}: {line}", line_num = line_idx + 1))
                        .map_err(Into::into),
                );
            }
        }

        None
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
/// - The ranges must not be perfectly adjacent (all ranges must either overlap
///   or have a non-coding region between them)
/// - A single region of overlap cannot involve more than 2 ranges
fn parse_exon_coords(coords: &str) -> std::io::Result<Vec<ExonCoords>> {
    /// The maximum amount of overlap allowed between exons.
    ///
    /// SARS-CoV-2 requires -1 exon-to-exon frameshift with other viruses
    /// reported up to -2.
    const MAX_DUPLICATED_OVERLAP_NT: usize = 2;

    let mut exon_ranges: Vec<ExonCoords> = Vec::new();
    let mut cds_start = 0;

    // Parses and pushes a range to exon_ranges. Iterators do not work since we
    // need to access the previous ExonCoords to validate order and overlap.
    for ref_range in parse_coords_inclusive::<Range<usize>>(coords) {
        let ref_range = ref_range?;

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

            // Prevent perfectly adjacent exons (there should either be overlap
            // or non-coding region)
            if last.ref_range.end == ref_range.start {
                return Err(std::io::Error::other(format!(
                    "Two exons are perfectly adjacent, and should therefore be represented as a single exon. Found {} then {}",
                    last.ref_range.display_inclusive(),
                    ref_range.display_inclusive(),
                )));
            }

            // Prevent overlapping exons that overlap by more than
            // MAX_DUPLICATED_OVERLAP_NT
            //
            // Exclusive index - inclusive index is valid length
            let overlap_nt = last.ref_range.end.saturating_sub(ref_range.start);
            if overlap_nt > MAX_DUPLICATED_OVERLAP_NT {
                return Err(std::io::Error::other(format!(
                    "Exon overlap exceeds {MAX_DUPLICATED_OVERLAP_NT} nt! Found {} then {}",
                    last.ref_range.display_inclusive(),
                    ref_range.display_inclusive(),
                )));
            }

            // Prevent a single region of overlap from involving more than 2
            // exons
            if let Some(second_to_last) = exon_ranges.iter().nth_back(1) {
                let overlap_nt = second_to_last.ref_range.end.saturating_sub(ref_range.start);
                if overlap_nt > 0 {
                    return Err(std::io::Error::other(format!(
                        "A single region of overlap cannot involve more than 2 exons within a given protein product. Found {}, {}, then {}, which all overlap",
                        second_to_last.ref_range.display_inclusive(),
                        last.ref_range.display_inclusive(),
                        ref_range.display_inclusive(),
                    )));
                }
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

trait BufReaderExt {
    fn starts_with(&mut self, needle: &[u8]) -> std::io::Result<bool>;
}

impl<R: Read> BufReaderExt for BufReader<R> {
    fn starts_with(&mut self, needle: &[u8]) -> std::io::Result<bool> {
        Ok(self.peek(needle.len())? == needle)
    }
}
