//! Codon position weight matrix loading.

use crate::data::keys::SpecKey;
use std::{
    cmp::Ordering,
    collections::{HashMap, HashSet, hash_map::Entry},
    fs::File,
    io::BufRead,
    path::Path,
    str::FromStr,
    sync::LazyLock,
};
use zoe::{
    data::{
        ByteMap,
        err::{ResultWithErrorContext, WithErrorContext},
    },
    iter_utils::ProcessResultsExt,
    prelude::{Nucleotides, NucleotidesView},
};

/// A map from spec keys to codon-position weights.
pub type CodonWeightMatrix = HashMap<SpecKey, CodonPositionWeights>;

/// A map counting the occurrence of different codons at each position in the
/// sequence family. The counts are also called _weights_.
///
/// This is stored using a [`HashMap`], which supports fast look-up although is
/// less cache-friendly. Since accesses into this data structure are considered
/// uncommon (to handle individual insertions or deletions), as compared to
/// looking up statistics for each position consecutively, this is sufficient.
///
/// This map assumes that all codons are in uppercase. Any `U` bases in the
/// codons are automatically converted to `T`.
#[derive(Clone, Debug, Default)]
pub struct CodonPositionWeights {
    map: HashMap<CodonKey, u32>,
}

impl CodonPositionWeights {
    /// Creates a new empty [`CodonPositionWeights`] map.
    pub fn new() -> Self {
        Self { map: HashMap::new() }
    }

    /// Returns whether the [`CodonPositionWeights`] map is empty.
    pub fn is_empty(&self) -> bool {
        self.map.is_empty()
    }

    /// Retrieves the count for a provided `codon` at the 1-based `position`. If
    /// it does not exist, 0 is returned.
    ///
    /// ## Validity
    ///
    /// The `codon` must contain unaligned, uppercase IUPAC bases. `T` must be
    /// used instead of `U`.
    pub fn get(&self, position: u32, codon: [u8; 3]) -> u32 {
        self.map.get(&CodonKey { position, codon }).copied().unwrap_or(0)
    }

    /// Inserts a 1-based position and codon into the map, returning the old
    /// count if present.
    ///
    /// ## Errors
    ///
    /// The position/codon pair must not already be present in the
    /// [`CodonPositionWeights`].
    ///
    /// ## Validity
    ///
    /// The `codon` must contain unaligned, uppercase IUPAC bases. `T` must be
    /// used instead of `U`.
    pub fn insert(&mut self, position: u32, codon: [u8; 3], count: u32) -> std::io::Result<()> {
        match self.map.entry(CodonKey { position, codon }) {
            Entry::Occupied(_) => Err(std::io::Error::other(format!(
                "The position/codon combination already exists in the specs.\n | Position: {position}\n | Codon: {codon}",
                codon = NucleotidesView::from(&codon)
            ))),
            Entry::Vacant(entry) => {
                entry.insert(count);
                Ok(())
            }
        }
    }

    /// Compares two counts of two codons at a specified 1-based position,
    /// returning `None` if both have count 0.
    ///
    /// ## Validity
    ///
    /// The `codon` must contain unaligned, uppercase IUPAC bases. `T` must be
    /// used instead of `U`.
    pub fn compare_codons(&self, left: [u8; 3], right: [u8; 3], position: u32) -> Option<Ordering> {
        let left_count = self.get(position, left);
        let right_count = self.get(position, right);
        if left_count == 0 && right_count == 0 {
            None
        } else {
            Some(left_count.cmp(&right_count))
        }
    }

    /// Compares two counts of the same codon at different 1-based positions,
    /// returning `None` if both have count 0.
    ///
    /// ## Validity
    ///
    /// The `codon` must contain unaligned, uppercase IUPAC bases. `T` must be
    /// used instead of `U`.
    pub fn compare_positions(&self, pos_left: u32, pos_right: u32, codon: [u8; 3]) -> Option<Ordering> {
        let left_count = self.get(pos_left, codon);
        let right_count = self.get(pos_right, codon);
        if left_count == 0 && right_count == 0 {
            None
        } else {
            Some(left_count.cmp(&right_count))
        }
    }
}

/// A key for codon position weights, combining position and the codon.
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug)]
struct CodonKey {
    /// The 1-based position of the codon within the coding sequence.
    position: u32,
    /// The uppercase codon, containing `T` instead of `U`.
    codon:    [u8; 3],
}

/// Default codon usage statistics for influenza A.
///
/// We convert codon usage statistics from
/// <http://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=465364&aa=1&style=GCG>
/// into ranks (times 10, in case we need to make easy adjustments), since that
/// is all that needed for our algorithm. Stop codons are set to the lowest
/// ranks to avoid premature stop codons.
#[rustfmt::skip]
pub(crate) static DEFAULT_CODON_STATS: LazyLock<HashMap<[u8; 3], u32>> = LazyLock::new(|| {
    HashMap::from([
        (*b"TAG",  10), (*b"TGA",  20), (*b"TAA",  30), (*b"CGT",  40),
        (*b"CGC",  50), (*b"TCG",  60), (*b"GCG",  70), (*b"CCG",  80),
        (*b"ACG",  90), (*b"CGG", 100), (*b"CGA", 110), (*b"CAC", 120),
        (*b"CCC", 130), (*b"TGT", 140), (*b"TTA", 150), (*b"GGC", 160),
        (*b"CAT", 170), (*b"GGT", 180), (*b"TCC", 190), (*b"GTC", 200),
        (*b"TGC", 210), (*b"CCT", 220), (*b"CTC", 230), (*b"GTA", 240),
        (*b"TAC", 250), (*b"GCC", 260), (*b"TCT", 270), (*b"ACC", 280),
        (*b"CTA", 290), (*b"GTT", 300), (*b"TAT", 310), (*b"AGT", 320),
        (*b"TTG", 330), (*b"CCA", 340), (*b"CTG", 350), (*b"AGC", 360),
        (*b"GCT", 370), (*b"TGG", 380), (*b"CTT", 390), (*b"TTT", 400),
        (*b"GGG", 410), (*b"ATC", 420), (*b"AGG", 430), (*b"ACT", 440),
        (*b"CAG", 450), (*b"TCA", 460), (*b"GTG", 470), (*b"TTC", 480),
        (*b"GAC", 490), (*b"AAG", 500), (*b"CAA", 510), (*b"AAC", 520),
        (*b"ATT", 530), (*b"GCA", 540), (*b"ATA", 550), (*b"GAT", 560),
        (*b"ACA", 570), (*b"GGA", 580), (*b"GAG", 590), (*b"AAT", 600),
        (*b"AGA", 610), (*b"AAA", 620), (*b"ATG", 630), (*b"GAA", 640),
    ])
});

/// Loads the codon position weights from a TSV file.
///
/// The file format uses pipe-delimited headers to start new sections:
///
/// ```text
/// #reference_id|protein
/// position<TAB>codon<TAB>count
/// ...
/// ```
///
/// If the file is empty, then the resulting [`CodonWeightMatrix`] will also be
/// empty.
///
/// ## Errors
///
/// All IO errors are propagated without path context. Headers (lines containing
/// `|`) must successfully parse. If any headers fail to parse, an error is
/// returned with context including the line and the expected format. Each TSV
/// row must parse successfully, and there should be no duplicate position/codon
/// pairs. Codons must only contains unaligned IUPAC bytes.
pub fn load_codon_weights(path: &Path) -> std::io::Result<CodonWeightMatrix> {
    let file = File::open(path)?;
    let reader = std::io::BufReader::new(file);
    CodonWeightMatrixParser::from_reader(reader).map(|out| out.0)
}

/// A parsed data from a single row of the `codon-position-weights.tsv` file for
/// the module.
struct TsvRow {
    /// The 1-based position of the codon within the coding sequence of the
    /// product.
    position:  u32,
    /// The uppercase IUPAC codon, with `U` substituted for `T`.
    codon:     [u8; 3],
    /// The count of the codon at the specified position.
    count:     u32,
    /// The original codon, before recoding. This is useful for error messages.
    raw_codon: [u8; 3],
}

impl FromStr for TsvRow {
    type Err = std::io::Error;

    /// Parses a `line` to a [`TsvRow`]. The string must be trimmed and
    /// non-empty, and should not be a header row.
    ///
    /// The codons are converted to uppercase, and `U` is substituted for `T`.
    ///
    /// ## Errors
    ///
    /// If any of the columns are missing or empty, then an error is returned.
    /// The position and count must successfully parse as `u32`, and the codon
    /// must be three characters. The position must also be nonzero.
    fn from_str(line: &str) -> Result<Self, Self::Err> {
        // End the iterator on empty fields, so that missing field errors appear
        let mut parts = line.split_ascii_whitespace().take_while(|s| !s.is_empty());

        // Get fields as string slices
        let Some(position) = parts.next() else {
            return Err(std::io::Error::other("Missing position field (first field)"));
        };
        let Some(codon) = parts.next() else {
            return Err(std::io::Error::other("Missing codon field (second field)"));
        };
        let Some(count) = parts.next() else {
            return Err(std::io::Error::other("Missing count field (third field)"));
        };

        let position = position
            .parse::<u32>()
            .with_context("Failed to parse position field (first field)")?;

        if position == 0 {
            return Err(std::io::Error::other("Found a position of 0, but position is 1-based"));
        }

        let raw_codon: [u8; 3] = codon.as_bytes().try_into().map_err(|_| {
            std::io::Error::other(format!(
                "The codon was expected to have 3 ASCII characters, but {} were found",
                codon.len()
            ))
        })?;

        let codon = normalize_codon(raw_codon)
            .with_context(format!("The codon {codon} is invalid", codon = Nucleotides::from(raw_codon)))?;

        let count = count
            .parse::<u32>()
            .with_context("Failed to parse count field (third field)")?;

        Ok(Self {
            position,
            codon,
            count,
            raw_codon,
        })
    }
}

/// Normalizes the codon for use with [`CodonWeightMatrix`].
///
/// This converts to uppercase IUPAC, including converting `U` to `T`.
///
/// ## Errors
///
/// An error is returned if a non-IUPAC character is encountered.
fn normalize_codon(mut codon: [u8; 3]) -> std::io::Result<[u8; 3]> {
    /// Filters non-IUPAC bytes, and uppercases. Converts `U` to `T`.
    const SANITIZE: ByteMap = ByteMap::all(0)
        .preserve(b"ACGTRYSWKMBDHVN")
        .map(b"acgtryswkmbdhvn", b"ACGTRYSWKMBDHVN")
        .map_to_one(b"uU", b'T');

    for base in &mut codon {
        let new_base = SANITIZE[*base];

        if new_base == 0 {
            return Err(std::io::Error::other(format!(
                "The character {base} is not permitted in a codon"
            )));
        } else {
            *base = new_base;
        }
    }

    Ok(codon)
}

/// A wrapper type around [`CodonWeightMatrix`] used for parsing.
struct CodonWeightMatrixParser(CodonWeightMatrix);

impl CodonWeightMatrixParser {
    fn new() -> Self {
        Self(HashMap::new())
    }

    /// Determines whether a line corresponds to a section header by determining
    /// whether it contains `|`. Leading `#` are stripped.
    ///
    /// ## Errors
    ///
    /// The header must contain at least two non-empty parts separated by `|`.
    fn parse_section_header(line: &str) -> std::io::Result<Option<SpecKey>> {
        line.contains('|')
            .then(|| {
                // Parse new section header
                let line = line.trim_start_matches('#').trim_ascii();
                let mut parts = line.split('|').map(str::trim_ascii).take_while(|str| !str.is_empty());

                let (Some(reference_id), Some(product_name)) = (parts.next(), parts.next()) else {
                    return Err(std::io::Error::other(format!(
                        "Invalid weight matrix header. Expected reference_id|protein, found {line}"
                    )));
                };

                Ok(SpecKey {
                    reference_id: reference_id.to_string(),
                    product_name: product_name.to_string(),
                })
            })
            .transpose()
    }

    /// Adds the data for the given section to `self`.
    ///
    /// ## Errors
    ///
    /// The given [`SpecKey`] must not already be in the hashmap. No sections
    /// can be non-empty, and no sections can contain duplicate codon-position
    /// pairs.
    fn process_section(&mut self, key: SpecKey, rows: impl Iterator<Item = TsvRow>) -> std::io::Result<()> {
        match self.0.entry(key) {
            Entry::Occupied(entry) => {
                let SpecKey {
                    reference_id,
                    product_name,
                } = entry.key();

                Err(std::io::Error::other(format!(
                    "A duplicate section for reference ID {reference_id} and product {product_name} was found"
                )))
            }
            Entry::Vacant(entry) => {
                let mut weights = CodonPositionWeights::new();
                let mut raw_codons = HashSet::new();

                for row in rows {
                    let insert_res = weights.insert(row.position, row.codon, row.count);

                    if let Err(e) = insert_res {
                        let normalized = row.codon;

                        let other_raw_codon = raw_codons.into_iter().find(|other_raw| {
                            normalize_codon(*other_raw).is_ok_and(|other_normalized| other_normalized == normalized)
                        });

                        let Some(other_raw_codon) = other_raw_codon else {
                            return Err(e);
                        };

                        if row.raw_codon != other_raw_codon {
                            return Err(e
                                .with_context(format!(
                                    "Both {} and {} were found in the specs",
                                    NucleotidesView::from(&other_raw_codon),
                                    NucleotidesView::from(&row.raw_codon)
                                ))
                                .into());
                        }

                        return Err(e);
                    }

                    raw_codons.insert(row.raw_codon);
                }

                if weights.is_empty() {
                    return Err(std::io::Error::other("No data was found under the header"));
                }

                entry.insert(weights);

                Ok(())
            }
        }
    }

    /// Parses the multi-sectioned TSV into `Self` from the given `reader`.
    ///
    /// ## Errors
    ///
    /// See [`parse_section_header`], [`TsvRow::from_str`], and
    /// [`process_section`]. Context including the line number and failing line
    /// is added for parsing errors.
    ///
    /// [`parse_section_header`]: CodonWeightMatrixParser::parse_section_header
    /// [`process_section`]: CodonWeightMatrixParser::process_section
    fn from_reader<R: BufRead>(reader: R) -> std::io::Result<Self> {
        let mut out = Self::new();

        reader.lines().process_results(|iter| {
            // Enumerate the lines for error messages, trim & filter empty, and
            // make it peekable for detecting headers
            let mut lines = iter
                .enumerate()
                .filter_map(|(line_idx, line)| {
                    let line = line.trim_ascii();
                    (!line.is_empty()).then_some((line_idx, line.to_owned()))
                })
                .peekable();

            // Empty file is ok
            let Some((line_idx, mut current_header_line)) = lines.next() else {
                return Ok(());
            };

            // Parse header for first line
            let Some(mut current_header) = Self::parse_section_header(&current_header_line)? else {
                return Err(std::io::Error::other(format!(
                    "Failed to parse section header at line {line_num}: {current_header_line}",
                    line_num = line_idx + 1
                )));
            };

            loop {
                // The next header and line, if found
                let mut next_header = None;

                // Get iterator of all regular lines under the current header,
                // aborting when a header is found and storing it to next_header
                let iter = std::iter::from_fn(|| {
                    let (line_idx, line) = lines.next()?;

                    match Self::parse_section_header(&line) {
                        Ok(Some(header)) => {
                            next_header = Some((header, line));
                            None
                        }
                        Ok(None) => Some(line.parse().with_context(format!(
                            "Failed to parse TSV data on line {line_num}: {line}",
                            line_num = line_idx + 1
                        ))),
                        Err(e) => Some(Err(e.with_context(format!(
                            "Failed to parse section header at line {line_num}: {line}",
                            line_num = line_idx + 1
                        )))),
                    }
                });

                // Process the header + lines
                iter.process_results(|iter| out.process_section(current_header, iter))?
                    .with_context(format!("Failed to parse the data for section: {current_header_line}"))?;

                match next_header {
                    Some((header, header_line)) => {
                        current_header = header;
                        current_header_line = header_line;
                    }
                    None => break,
                }
            }

            Ok(())
        })??;

        Ok(out)
    }
}
