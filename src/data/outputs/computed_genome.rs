use crate::data::ranges::InsertionIdx;
use zoe::prelude::*;

/// Pre-computed genome data for `.gen` output.
///
/// This struct holds the expensive-to-compute data for genome output.
#[derive(Debug)]
pub struct PrecomputedGenomeData {
    /// SHA1 hash of genome sequence, or `None` if no DNA data remained after
    /// filtering.
    pub genome_id:     Option<String>,
    /// Genome length (ungapped)
    pub genome_length: usize,
    /// Whether any insertion exists
    pub has_insertion: bool,
    /// Genome sequence (without deletions, includes insertions)
    pub genome_seq:    Nucleotides,
    /// Genome alignment (with `-` for deletions, no insertions)
    pub genome_aln:    Nucleotides,
    /// Genome-level insertions
    pub insertions:    Vec<ComputedGenomeInsertion>,

    pub trailing_ref_unaligned: usize,
}

/// Genome-level insertion for `.gen.ins` output.
#[derive(Debug)]
pub struct ComputedGenomeInsertion {
    /// Upstream nucleotide position (1-based)
    pub upstream_nt_pos:      usize,
    /// Inserted nucleotides
    pub inserted_nucleotides: Nucleotides,
}

impl ComputedGenomeInsertion {
    pub fn new(nt_insertion_idx: InsertionIdx, slice: &[u8]) -> Self {
        ComputedGenomeInsertion {
            upstream_nt_pos:      nt_insertion_idx.left_pos(),
            inserted_nucleotides: Nucleotides::from(slice),
        }
    }
}
