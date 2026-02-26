use std::cmp::Ordering;

use crate::data::{
    exons::Exons,
    products::Product,
    ranges::StateRange,
    weights::{CodonPositionWeights, DEFAULT_CODON_STATS},
};

/// Protein product specification.
#[derive(Debug)]
pub(crate) struct ProductSpec {
    /// Protein or peptide name
    pub(crate) name:          String,
    /// Exon Coordinates and translation rules
    pub(crate) exons:         Exons,
    /// The protein product specific codon positive weight matrix
    pub(crate) codon_weights: Option<CodonPositionWeights>,
}

impl ProductSpec {
    pub(crate) fn make_product_ranges<'a>(&'a self, state_ranges: &[StateRange]) -> Product<'a> {
        // TODO: Is this a good enough capacity? We could end up exceeding it.
        let mut product_ranges = Vec::with_capacity(state_ranges.len());

        for state in state_ranges {
            for exon in &self.exons.coords {
                // A state can span multiple exons: e.g., long match for full contig
                // An exon can span multiple states: e.g., a match with an indel
                if let Some(cds_state) = state.intersect_exon(exon) {
                    product_ranges.push(cds_state);
                }
            }
        }

        Product {
            product_ranges,
            product_spec: self,
            stop_extension_query_range: None,
        }
    }

    /// Compares the counts of two codons at the specified 1-based position,
    /// returning true if `left >= right`.
    ///
    /// This uses position-specific weights if available, falling back to
    /// [`DEFAULT_CODON_STATS`] when both position-specific counts are zero.
    ///
    /// ## Validity
    ///
    /// The `left` and `right` codons must contain unaligned, uppercase IUPAC
    /// bases.
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
    /// The `codon` must contain unaligned, uppercase IUPAC bases.
    pub(crate) fn compare_codon_positions(&self, pos_left: u32, pos_right: u32, codon: [u8; 3]) -> Option<Ordering> {
        self.codon_weights
            .as_ref()
            .and_then(|w| w.compare_positions(pos_left, pos_right, codon))
    }
}
