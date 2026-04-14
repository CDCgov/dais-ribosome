//! The [`ProductSpec`] specification struct, containing references to a protein
//! product's name, exons, and codon position weights.

use crate::{
    data::{
        exons::Exons,
        ranges::{CdsStateRange, StateRange},
        weights::{CodonPositionWeights, DEFAULT_CODON_STATS},
    },
    error,
    outputs::Product,
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
    /// Intersects the ranges for an alignment ([`StateRange`]) with the ranges
    /// for the exons ([`Exons`]) to form the ranges in the product.
    ///
    /// The `stop_extension_query_range` field is initialized to `None`, and
    /// must be updated later.
    ///
    /// ## Validity
    ///
    /// The `state_ranges` must contain ordered non-overlapping ranges that
    /// fully partition the aligned query and reference ranges (if any part of
    /// the sequences was not locally aligned, this will not be included). It
    /// also must begin and end with [`StateRange::M`] (and hence be non-empty).
    pub(crate) fn make_product_ranges(&self, state_ranges: &[StateRange]) -> Product<'_> {
        let product_ranges = self.intersect(state_ranges);

        let leading_cds_unaligned = self.compute_leading_cds_unaligned(&product_ranges);
        let trailing_cds_unaligned = self.compute_trailing_cds_unaligned(&product_ranges);

        Product {
            product_spec: self,
            product_ranges,
            leading_cds_unaligned,
            trailing_cds_unaligned,
            stop_extension_query_range: None,
        }
    }

    /// A helper function for [`make_product_ranges`] which computes the
    /// `leading_cds_unaligned` field.
    ///
    /// [`make_product_ranges`]: ProductSpec::make_product_ranges
    fn compute_leading_cds_unaligned(&self, product_ranges: &[CdsStateRange]) -> usize {
        let Some(first) = product_ranges.first() else {
            // No intersection between query and exons, so put all the CDS
            // length as leading unaligned
            return self.exons.cds_len();
        };

        match first {
            CdsStateRange::M(m) => m.cds_range.start,
            CdsStateRange::D(d) => d.cds_range.start,
            CdsStateRange::I(i) => {
                error!("product_ranges cannot begin with an insertion!");
                i.cds_index.right()
            }
        }
    }

    /// A helper function for [`make_product_ranges`] which computes the
    /// `product_ranges` field.
    ///
    /// [`make_product_ranges`]: ProductSpec::make_product_ranges
    fn intersect(&self, state_ranges: &[StateRange]) -> Vec<CdsStateRange> {
        // TODO: Is this a good enough capacity? We could end up exceeding it.
        let mut product_ranges = Vec::with_capacity(state_ranges.len());

        // We want product_ranges ordered by coding sequence coordinates so that
        // product-specific edits to the alignment can be made (e.g., index
        // shifting). This is different from ordering by query coordinates,
        // which may have a differing order if any exons overlap. As such, the
        // outer loop is over the exons (coding sequence coordinates) and the
        // inner loop is over the state ranges (query coordinates).
        for exon in &self.exons.coords {
            for state in state_ranges {
                // A state can span multiple exons, such as a long match for a
                // full contig. An exon can span multiple states, such as a
                // match with an indel
                if let Some(cds_state) = state.intersect_exon(exon) {
                    product_ranges.push(cds_state);
                }
            }
        }

        product_ranges
    }

    /// A helper function for [`make_product_ranges`] which computes the
    /// `trailing_cds_unaligned` field.
    ///
    /// [`make_product_ranges`]: ProductSpec::make_product_ranges
    fn compute_trailing_cds_unaligned(&self, product_ranges: &[CdsStateRange]) -> usize {
        let Some(last) = product_ranges.last() else {
            // No intersection between query and exons, and we counted all the
            // CDS length as leading unaligned
            return 0;
        };

        match last {
            CdsStateRange::M(m) => m.cds_range.end,
            CdsStateRange::D(d) => d.cds_range.end,
            CdsStateRange::I(i) => {
                error!("product_ranges cannot begin with an insertion!");
                i.cds_index.right()
            }
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
