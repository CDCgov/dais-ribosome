//! Logic for intersecting alignment ranges ([`StateRange`]) with the exons
//! ([`Exons`]) in order to form a [`Product`].

use crate::{
    config::ProductSpec,
    data::exons::{ExonCoords, Exons},
    outputs::Product,
    ranges::{
        CdsDeletionRange, CdsInsertionRange, CdsMatchRange, CdsStateRange, DeletionRange, InsertionIdx, InsertionRange,
        MatchRange, RangeExt, StateRange,
    },
};

/// Intersects the ranges for an alignment ([`StateRange`]) with the ranges for
/// the exons ([`Exons`]) to form the ranges in the product.
///
/// The `stop_extension_query_range` field is initialized to `None`, and must be
/// updated later.
///
/// ## Validity
///
/// The `state_ranges` must contain ordered non-overlapping ranges that fully
/// partition the aligned query and reference ranges (if any part of the
/// sequences was not locally aligned, this will not be included). None of the
/// ranges can be empty. It also must begin and end with [`StateRange::M`].
pub(crate) fn form_product<'a>(state_ranges: &[StateRange], product_spec: &'a ProductSpec) -> Product<'a> {
    let product_ranges = intersect_with_exons(state_ranges, &product_spec.exons);

    let leading_cds_unaligned = match product_ranges.first() {
        Some(CdsStateRange::M(m)) => m.cds_range.start,
        Some(CdsStateRange::D(d)) => d.cds_range.start,
        Some(CdsStateRange::I(i)) => {
            // Unreachable in the current algorithm but may change in the future
            i.cds_index.right()
        }

        // We put all of the unaligned bases in trailing_cds_unaligned
        None => 0,
    };

    let trailing_cds_unaligned = {
        let end = match product_ranges.last() {
            Some(CdsStateRange::M(m)) => m.cds_range.end,
            Some(CdsStateRange::D(d)) => d.cds_range.end,
            Some(CdsStateRange::I(i)) => {
                // Unreachable in the current algorithm but may change in the
                // future
                i.cds_index.right()
            }

            // If product_ranges is empty, then the aligned-against region ends
            // at 0 (resulting in trailing_cds_unaligned being cds_len)
            None => 0,
        };

        product_spec.exons.cds_len() - end
    };

    Product {
        product_spec,
        product_ranges,
        leading_cds_unaligned,
        trailing_cds_unaligned,
        stop_extension_query_range: None,
    }
}

/// A helper function for [`form_product`] which computes the `product_ranges`
/// field.
///
/// The output `Vec` may be empty (in which case there was no intersection of
/// the state ranges with any exons). It may begin and end with match states or
/// delete states, but it will never begin or end with an insert state.
///
/// Any returned [`CdsStateRange`] entries will have a non-zero length.
///
/// ## Validity
///
/// Same validity requirements as [`form_product`].
fn intersect_with_exons(state_ranges: &[StateRange], exons: &Exons) -> Vec<CdsStateRange> {
    // This capacity may be exceeded, but is a good initial choice
    let mut product_ranges = Vec::with_capacity(state_ranges.len());

    // We want product_ranges ordered by coding sequence coordinates so that
    // product-specific edits to the alignment can be made (e.g., frame
    // shifting). This is different from ordering by query coordinates, which
    // may have a differing order if any exons overlap. As such, the outer loop
    // is over the exons (coding sequence coordinates) and the inner loop is
    // over the state ranges (query coordinates).
    for exon in &exons.coords {
        for state in state_ranges {
            // A state can span multiple exons, such as a long match for a
            // full contig. An exon can span multiple states, such as a
            // match with an indel
            if let Some(cds_state) = state.intersect_exon(exon) {
                product_ranges.push(cds_state);
            }
        }
    }

    // Validity: the output range will not begin or end with an insertion. An
    // intersected insertion implies that the insertion is strictly contained in
    // an exon, but therefore there must exist flanking ranges in `state_ranges`
    // that also intersect the exon. This is because `state_ranges` forms an
    // ordered partition of the aligned query and reference ranges, and cannot
    // begin or end with an insertion.
    product_ranges
}

impl MatchRange {
    /// Intersects a [`MatchRange`] with an exon in reference coordinates,
    /// returning the resulting query coordinates and coding sequence
    /// coordinates.
    ///
    /// If `Some`, the [`CdsMatchRange`] will have a non-zero length.
    fn intersect_exon(&self, exon: &ExonCoords) -> Option<CdsMatchRange> {
        // Intersect the two ranges in reference coordinates
        self.ref_range.intersect(&exon.ref_range).map(|intersect_ref_range| {
            let self_shrinkage = self.ref_range.compute_shrinkage(&intersect_ref_range);
            let intersect_query_range = self.query_range.shrink(self_shrinkage);

            // self.query_range and self.ref_range have the same length, so the
            // same should be true for the intersected versions
            debug_assert_eq!(intersect_query_range.len(), intersect_ref_range.len());

            let exon_shrinkage = exon.ref_range.compute_shrinkage(&intersect_ref_range);
            let intersect_cds_range = exon.cds_range.shrink(exon_shrinkage);

            // exon.cds_range and exon.ref_range have the same length, so the
            // same should be true for the intersected versions
            debug_assert_eq!(intersect_cds_range.len(), intersect_ref_range.len());

            // Validity: Per above, these are the same length, equal to the
            // intersect_ref_range.len(). This is non-zero due to intersect
            // guarantees
            CdsMatchRange {
                query_range: intersect_query_range,
                cds_range:   intersect_cds_range,
            }
        })
    }
}

impl DeletionRange {
    /// Intersects a [`DeletionRange`] with an exon in reference coordinates,
    /// returning the resulting coding sequence coordinates.
    ///
    /// If `Some`, the [`CdsMatchRange`] will have a non-zero length.
    fn intersect_exon(&self, exon: &ExonCoords) -> Option<CdsDeletionRange> {
        self.ref_range.intersect(&exon.ref_range).map(|intersect_ref_range| {
            let exon_shrinkage = exon.ref_range.compute_shrinkage(&intersect_ref_range);
            let intersect_cds_range = exon.cds_range.shrink(exon_shrinkage);

            // Validity: exon.ref_range and exon.cds_range are the same length.
            // So, intersect_cds_range is the same length as
            // intersect_ref_range, which is non-empty due to intersect
            // guarantees
            CdsDeletionRange {
                cds_range: intersect_cds_range,
            }
        })
    }
}

impl InsertionRange {
    /// If the insertion range and exon strictly intersect (the insertion
    /// appears in the middle of the exon), then compute the
    /// [`CdsInsertionRange`] of the intersection.
    ///
    /// ## Validity
    ///
    /// The length of the insertion must be non-zero in order to ensure the
    /// output has a non-zero length.
    fn intersect_exon(&self, exon: &ExonCoords) -> Option<CdsInsertionRange> {
        exon.ref_range.contains_ins(self.ref_index).then(|| {
            // Intuitively, cds_range.start + (ref_index.right-ref_range.start)
            // where the second term is the offset of the insertion within the
            // reference range. However, that offset may be positive or
            // negative, so we need to do additions before substractions to
            // prevent underflow
            let cds_index =
                InsertionIdx::from_right_idx(exon.cds_range.start + self.ref_index.right() - exon.ref_range.start);

            CdsInsertionRange {
                cds_index,
                query_range: self.query_range.clone(),
            }
        })
    }
}

impl StateRange {
    /// Intersects a [`StateRange`] with an exon in reference coordinates,
    /// returning the resulting query coordinates and coding sequence
    /// coordinates.
    ///
    /// If the return value is `Some`, the [`CdsStateRange`] will have a
    /// non-zero length.
    ///
    /// An insertion is only considered to intersect an exon if it occurs
    /// strictly within the exon (not on the boundaries).
    ///
    /// ## Validity
    ///
    /// The contained ranges in `self` must be non-empty.
    fn intersect_exon(&self, exon: &ExonCoords) -> Option<CdsStateRange> {
        match self {
            Self::M(m) => m.intersect_exon(exon).map(CdsStateRange::M),
            Self::D(d) => d.intersect_exon(exon).map(CdsStateRange::D),
            Self::I(i) => i.intersect_exon(exon).map(CdsStateRange::I),
        }
    }
}
