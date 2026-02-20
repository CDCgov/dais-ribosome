//! Codon position weight matrix loading.

use crate::data::{
    TSVReader,
    keys::{CodonKey, SpecKey},
};
use std::{
    collections::HashMap,
    io::{BufRead, BufReader},
    path::Path,
    sync::LazyLock,
};

/// Default codon usage statistics for Influenza A
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

/// Map from codon key (position + codon) to weight/count.
pub type CodonPositionWeights = HashMap<CodonKey, u32>;

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
    let file = TSVReader::open_nonempty_file(path)?;
    let reader = BufReader::new(file);
    let mut all_matrices = HashMap::new();
    let mut current_weights = HashMap::new();
    let mut current_key: Option<SpecKey> = None;

    for line in reader.lines() {
        let line = line?;
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
        }

        // Parse weight entry: position<TAB>codon<TAB>count
        let mut parts = line.split_ascii_whitespace();
        if let (Some(position_str), Some(codon_str), Some(count_str)) = (parts.next(), parts.next(), parts.next())
            && let (Ok(position), Ok(count), Some(Ok(codon))) = (
                position_str.trim_ascii().parse(),
                count_str.trim_ascii().parse::<u32>(),
                codon_str.trim_ascii().as_bytes().get(..3).map(TryInto::try_into),
            )
        {
            let mut key = CodonKey { position, codon };
            key.codon.make_ascii_uppercase();
            current_weights.insert(key, count);
        }
    }

    // Save final section
    if let Some(key) = current_key
        && !current_weights.is_empty()
    {
        all_matrices.insert(key, current_weights);
    }

    Ok(all_matrices)
}
