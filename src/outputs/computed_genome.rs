use crate::data::ranges::InsertionIdx;
use zoe::prelude::Nucleotides;

/// A computed genome, with materialized nucleotide and amino acid sequences.
#[derive(Debug)]
pub struct ComputedGenome {
    /// The SHA1 hash of the cleaned genome sequence, or `None` if no DNA data
    /// remained after filtering.
    pub genome_id:              Option<String>,
    /// The length of the genome's unaligned nucleotide sequence.
    ///
    /// This is equivalent to `genome_seq.len()`.
    pub genome_length:          usize,
    /// Whether any insertion exists in the genome.
    pub has_insertion:          bool,
    /// The unaligned nucleotide sequence for the genome (with insertions but no
    /// deletions).
    ///
    /// This will only contain unaligned uppercase IUPAC. Both `U` and `T` are
    /// allowed.
    pub genome_seq:             Nucleotides,
    /// The aligned nucleotide sequence for the genome (with `-` for deletions
    /// but no insertions).
    ///
    /// This will only contain uppercase IUPAC, padding `.`, and gaps `-`. Both
    /// `U` and `T` are allowed.
    pub genome_aln:             Nucleotides,
    /// The computed insertions within the genome.
    ///
    /// Unlike [`ComputedProduct::insertions`], no insertions are filtered.
    ///
    /// [`ComputedProduct::insertions`]:
    ///     crate::outputs::ComputedProduct::insertions
    pub insertions:             Vec<ComputedGenomeInsertion>,
    /// The number of bases in the reference sequence that were not aligned
    /// against at the end.
    ///
    /// See [`GenomeAndProductStates::trailing_ref_unaligned`].
    ///
    /// [`GenomeAndProductStates::trailing_ref_unaligned`]:
    ///     crate::outputs::GenomeAndProductStates::trailing_ref_unaligned
    pub trailing_ref_unaligned: usize,
}

/// Genome-level insertion for `.gen.ins` output.
#[derive(Debug)]
pub struct ComputedGenomeInsertion {
    /// The upstream nucleotide position (1-based), which is the position
    /// _after_ which the insertion occurs.
    pub upstream_nt_pos: usize,
    /// The inserted nucleotides.
    ///
    /// This will only contain unaligned uppercase IUPAC. Both `U` and `T` are
    /// allowed.
    pub inserted_nt:     Nucleotides,
}

impl ComputedGenomeInsertion {
    /// Creates a [`ComputedGenomeInsertion`] from raw insertion data and the
    /// insertion index.
    ///
    /// ## Validity
    ///
    /// The `inserted` bases must contain unaligned uppercase IUPAC. Both `U`
    /// and `T` are allowed. This is true for any slice of
    /// [`QueryRecord::nucleotides`].
    ///
    /// [`QueryRecord::nucleotides`]: crate::data::QueryRecord::nucleotides
    pub(crate) fn new(nt_insertion_idx: InsertionIdx, inserted: &[u8]) -> Self {
        ComputedGenomeInsertion {
            upstream_nt_pos: nt_insertion_idx.left_pos(),
            inserted_nt:     Nucleotides::from(inserted),
        }
    }
}
