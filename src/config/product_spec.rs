//! The [`ProductSpec`] specification struct, containing references to a protein
//! product's name, exons, and codon position weights.

use crate::data::{
    exons::Exons,
    weights::{CodonPositionWeights, DEFAULT_CODON_STATS},
};
use std::cmp::Ordering;

/// The specifications for a single protein product (e.g., `HA`, `HA-signal`).
#[derive(Clone, Debug)]
pub(crate) struct ProductSpec {
    /// The protein/peptide product name
    pub(crate) name:          String,
    /// The exon coordinates, as well as translation rules specific to the exons
    pub(crate) exons:         Exons,
    /// The codon position weight matrix for the protein product
    pub(crate) codon_weights: Option<CodonPositionWeights>,
}

impl ProductSpec {
    /// Compares the counts of two codons at the specified 1-based position,
    /// returning true if `left >= right`.
    ///
    /// This uses position-specific weights if available, falling back to
    /// [`DEFAULT_CODON_STATS`] when both position-specific counts are zero.
    ///
    /// ## Validity
    ///
    /// The `left` and `right` codons must contain unaligned, uppercase IUPAC
    /// bases. `T` must be used instead of `U`.
    pub(crate) fn codon_left_ge_right(&self, left: [u8; 3], right: [u8; 3], codon_position: u32) -> bool {
        // Validity: both codons are in uppercase
        if let Some(w) = &self.codon_weights
            && let Some(cmp) = w.compare_codons(left, right, codon_position)
        {
            cmp.is_ge()
        } else {
            let left_count = DEFAULT_CODON_STATS.get(&left).copied().unwrap_or(0);
            let right_count = DEFAULT_CODON_STATS.get(&right).copied().unwrap_or(0);
            left_count >= right_count
        }
    }

    // TODO: Use this!!

    /// Compares the likelihood of two codons at the specified 1-based position,
    /// returning an ordering based on the observed counts.
    ///
    /// If both observed counts are 0, or if there are no position-specific
    /// weights are available, [`DEFAULT_CODON_STATS`] is used.
    ///
    /// ## Validity
    ///
    /// The `left` and `right` codons must contain unaligned, uppercase IUPAC
    /// bases.
    #[allow(dead_code)]
    pub(crate) fn compare_codons(&self, left: [u8; 3], right: [u8; 3], position: u32) -> Ordering {
        // Validity: both codons are in uppercase
        self.codon_weights
            .as_ref()
            .and_then(|w| w.compare_codons(left, right, position))
            .unwrap_or_else(|| {
                let left_count = DEFAULT_CODON_STATS.get(&left).copied().unwrap_or(0);
                let right_count = DEFAULT_CODON_STATS.get(&right).copied().unwrap_or(0);
                left_count.cmp(&right_count)
            })
    }

    /// Compares the likelihood of a codon appearing at the 1-based positions
    /// `pos_left` and `pos_right`, returning an ordering based on the observed
    /// counts.
    ///
    /// If both observed counts are 0, or if there are no observed counts
    /// available, `None` is returned.
    ///
    /// ## Validity
    ///
    /// The `codon` must contain unaligned, uppercase IUPAC bases. `T` must be
    /// used instead of `U`.
    pub(crate) fn compare_codon_positions(&self, pos_left: u32, pos_right: u32, codon: [u8; 3]) -> Option<Ordering> {
        self.codon_weights
            .as_ref()
            .and_then(|w| w.compare_positions(pos_left, pos_right, codon))
    }
}
