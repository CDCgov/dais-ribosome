use std::cmp::Ordering;

use crate::data::{
    exons::Exons,
    products::Product,
    ranges::{CdsStateRange, StateRange},
    weights::{CodonPositionWeights, DEFAULT_CODON_STATS},
};

/// The specifications for a single protein product (e.g., `HA`, `HA-signal`).
#[derive(Copy, Clone, Debug)]
pub(crate) struct ProductSpec<'a> {
    /// The protein/peptide product name
    pub(crate) name:          &'a str,
    /// The exon coordinates, as well as translation rules specific to the exons
    pub(crate) exons:         &'a Exons,
    /// The codon position weight matrix for the protein product
    pub(crate) codon_weights: Option<&'a CodonPositionWeights>,
}

impl<'a> ProductSpec<'a> {
    /// Intersects the ranges for an alignment ([`StateRange`]) with the ranges
    /// for the exons ([`Exons`]) to form the ranges in the product.
    ///
    /// The `stop_extension_query_range` field is initialized to `None`, and
    /// must be updated later.
    ///
    /// ## Validity
    ///
    /// The `state_ranges` must contain ordered non-overlapping ranges that
    /// fully partition the query and reference ranges included in the
    /// alignment. It also must begin and end with [`StateRange::M`].
    pub(crate) fn make_product_ranges(self, state_ranges: &[StateRange]) -> Product<'a> {
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

        // Validity: the states and the exons are both ordered and
        // non-overlapping, so the above loop will maintain this property for
        // product_ranges. The first product range is not an insertion, since
        // this would imply that the range immediately before the insertion
        // (which exists since state_ranges starts with a match) should also
        // intersect the exon. Similarly for the end of the product_ranges.
        // Since state_ranges fully partition the reference sequence,
        // product_ranges will fully partition the exons, excluding exons at the
        // beginning/end which are not aligned against (or partially aligned
        // against)

        // TODO: Add validity that at least one match or delete must be present,
        // so that unwrap_or can be removed

        let leading_cds_unaligned = product_ranges
            .iter()
            .find_map(|s| match s {
                CdsStateRange::M(m) => Some(m.cds_range.start),
                CdsStateRange::D(d) => Some(d.cds_range.start),
                _ => None,
            })
            .unwrap_or(0);

        let trailing_cds_unaligned = self.exons.total_cds_length
            - product_ranges
                .iter()
                .rev()
                .find_map(|s| match s {
                    CdsStateRange::M(m) => Some(m.cds_range.end),
                    CdsStateRange::D(d) => Some(d.cds_range.end),
                    _ => None,
                })
                .unwrap_or(0);

        Product {
            product_spec: self,
            product_ranges,
            leading_cds_unaligned,
            trailing_cds_unaligned,
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
