//! Codon position weight matrix loading.

use crate::data::keys::{CodonKey, SpecKey};
use std::{
    cmp::Ordering,
    collections::HashMap,
    fs::File,
    io::{BufRead, BufReader},
    path::Path,
    str::FromStr,
    sync::LazyLock,
};
use zoe::data::err::ResultWithErrorContext;

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

/// A map counting the occurrence of different codons at each position in the
/// sequence family. The counts are also called _weights_.
///
/// This is stored using a [`HashMap`], which supports fast look-up although is
/// less cache-friendly. Since accesses into this data structure are considered
/// uncommon (to handle individual insertions or deletions), as compared to
/// looking up statistics for each position consecutively, this is sufficient.
///
/// This map assumes that all codons are in uppercase.
#[derive(Debug, Default)]
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
    /// The `codon` should be in uppercase.
    pub fn get(&self, position: u32, codon: [u8; 3]) -> u32 {
        self.map.get(&CodonKey { position, codon }).copied().unwrap_or(0)
    }

    /// Inserts a 1-based position and codon into the map, returning the old
    /// count if present.
    ///
    /// ## Validity
    ///
    /// The `codon` should be in uppercase.
    pub fn insert(&mut self, position: u32, codon: [u8; 3], count: u32) -> Option<u32> {
        self.map.insert(CodonKey { position, codon }, count)
    }

    /// Compares two counts of two codons at a specified 1-based position,
    /// returning `None` if both have count 0.
    ///
    /// ## Validity
    ///
    /// Both codons should be in uppercase.
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
    /// The codon should be in uppercase.
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

/// Map from spec key to codon position weights.
pub type CodonWeightMatrix = HashMap<SpecKey, CodonPositionWeights>;

/// Load codon position weight matrices from a TSV file.
///
/// The file format uses pipe-delimited headers to start new sections:
///
/// ```text
/// reference_id|protein
/// position<TAB>codon<TAB>count
/// ...
/// ```
///
/// ## Errors
///
/// Returns an error if:
///
/// - The file cannot be read
/// - A section header is malformed
pub fn load_codon_weights(path: &Path) -> Result<CodonWeightMatrix, std::io::Error> {
    // TODO: Check non-empty
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut all_matrices = HashMap::new();
    let mut current_weights = CodonPositionWeights::new();
    let mut current_key: Option<SpecKey> = None;

    for (line_idx, line) in reader.lines().enumerate() {
        let line = line?;
        let line = line.trim_ascii();

        if line.contains('|') {
            // Save previous section if it exists
            if let Some(key) = current_key.take()
                && !current_weights.is_empty()
            {
                all_matrices.insert(key, std::mem::take(&mut current_weights));
            }

            // Parse new section header
            let mut parts = line.trim_start_matches('#').trim_ascii().split('|');

            let (Some(reference_id), Some(protein_product)) = (parts.next(), parts.next()) else {
                return Err(std::io::Error::other("Issue parsing weight matrix header"));
            };

            current_key = Some(SpecKey::new(reference_id, protein_product));

            continue;
        } else if line.is_empty() {
            continue;
        }

        // Parse weight entry: position<TAB>codon<TAB>count
        let row = line
            .parse::<TsvRow>()
            .with_context(format!("Failed to parse line {line_num}: {line}", line_num = line_idx + 1))?;

        // Insert into the current weight map
        current_weights.insert(row.position, row.codon, row.count);
    }

    // Save final section
    if let Some(key) = current_key
        && !current_weights.is_empty()
    {
        all_matrices.insert(key, current_weights);
    }

    Ok(all_matrices)
}

struct TsvRow {
    position: u32,
    codon:    [u8; 3],
    count:    u32,
}

// TODO: We should validate that the code is only ACTG, or specify what
// behavior is otherwise. For example, RTC codon appears in
// flu-codon-position-weights.tsv, which will never match.

impl FromStr for TsvRow {
    type Err = std::io::Error;

    /// Parses a string `line` to a [`TsvRow`]. The string must be trimmed and
    /// non-empty, and should not be a header row.
    ///
    /// The codons are converted to uppercase.
    fn from_str(line: &str) -> Result<Self, Self::Err> {
        // End the iterator on empty fields, so that missing field errors appear
        let mut parts = line.split('\t').map(str::trim_ascii).take_while(|s| !s.is_empty());

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

        let mut codon: [u8; 3] = codon.as_bytes().try_into().map_err(|_| {
            std::io::Error::other(format!(
                "The codon was expected to have 3 ASCII characters, but {} were found",
                codon.len()
            ))
        })?;

        codon.make_ascii_uppercase();

        let count = count
            .parse::<u32>()
            .with_context("Failed to parse count field (third field)")?;

        Ok(Self { position, codon, count })
    }
}
