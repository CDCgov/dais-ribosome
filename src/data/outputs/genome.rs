use crate::{
    config::toml::Formatting,
    data::{Nullable, ranges::DeletionRange},
};
use std::fmt::{self, Display};
use zoe::prelude::*;

/// Genome-level insertion for `.gen.ins` output.
#[derive(Debug)]
pub struct ComputedGenomeInsertion {
    /// Upstream nucleotide position (1-based)
    pub upstream_nt:          usize,
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
}

/// A single row for `.gen` output.
pub struct GenRow<'a> {
    pub id:                &'a str,
    pub ctype:             &'a str,
    pub ref_id:            &'a str,
    pub genome:            &'a PrecomputedGenomeData,
    pub(crate) ref_len:    usize,
    pub(crate) formatting: &'a Formatting,
}

impl Display for GenRow<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let g = self.genome;
        let gen_rpad = if self.formatting.right_pad_gen {
            let trailing = self.ref_len.saturating_sub(g.genome_aln.len());
            ".".repeat(trailing)
        } else {
            String::new()
        };
        write!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}{}",
            self.id,
            self.ctype,
            self.ref_id,
            Nullable(&g.genome_id),
            g.genome_length,
            g.has_insertion,
            g.genome_seq,
            g.genome_aln,
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
            self.id, self.ctype, self.ref_id, ins.upstream_nt, ins.inserted_nucleotides,
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
