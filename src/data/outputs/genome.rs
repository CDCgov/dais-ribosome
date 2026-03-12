use crate::{
    config::Formatting,
    data::{Nullable, ranges::DeletionRange},
};
use std::fmt::{self, Display};
use zoe::prelude::*;

/// Genome-level insertion for `.gen.ins` output.
#[derive(Debug)]
pub struct ComputedGenomeInsertion {
    /// Upstream nucleotide position (1-based)
    pub upstream_nt_pos:      usize,
    /// Inserted nucleotides
    pub inserted_nucleotides: Nucleotides,
}

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

/// A single row for `.gen` output.
pub struct GenRow<'a> {
    pub id:                &'a str,
    pub ctype:             &'a str,
    pub ref_id:            &'a str,
    pub genome:            &'a PrecomputedGenomeData,
    pub(crate) formatting: &'a Formatting,
}

impl Display for GenRow<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let gen_rpad = if self.formatting.right_pad_gen {
            ".".repeat(self.genome.trailing_ref_unaligned)
        } else {
            String::new()
        };
        write!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}{}",
            self.id,
            self.ctype,
            self.ref_id,
            Nullable(&self.genome.genome_id),
            self.genome.genome_length,
            self.genome.has_insertion,
            self.genome.genome_seq,
            self.genome.genome_aln,
            gen_rpad,
        )
    }
}

/// A single row for `.gen.ins` output.
pub struct GenInsRow<'a> {
    pub id:        &'a str,
    pub ctype:     &'a str,
    pub ref_id:    &'a str,
    pub insertion: &'a ComputedGenomeInsertion,
}

impl Display for GenInsRow<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let ins = self.insertion;
        write!(
            f,
            "{}\t{}\t{}\t{}\t{}",
            self.id, self.ctype, self.ref_id, ins.upstream_nt_pos, ins.inserted_nucleotides,
        )
    }
}

/// A single row for `.gen.del` output.
pub struct GenDelRow<'a> {
    pub id:              &'a str,
    pub ctype:           &'a str,
    pub ref_id:          &'a str,
    pub(crate) deletion: &'a DeletionRange,
}

impl Display for GenDelRow<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let r = &self.deletion.ref_range;
        write!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{}",
            self.id,
            self.ctype,
            self.ref_id,
            r.start + 1,
            r.end,
            r.len(),
        )
    }
}
