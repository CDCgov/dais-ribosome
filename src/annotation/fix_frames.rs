//! Frame-fixing logic for indels.
//!
//! This module contains the functionality for fixing the frames of indels by
//! repositioning insertions and deletions to in-frame boundaries, as well as
//! merging adjacent deletions as long as they are not separated by noncoding
//! regions.
//!
//! This procedure does not physically reorder the bases in the query, but
//! rather reassigns which ones are matches and which ones are indels. For
//! example, when an insertion is shifted, it is the insertion operations in the
//! CIGAR string that are shifted, so that the insertion now contains different
//! bases.
//!
//! As such, this procedure may reduce the optimality of the Smith-Waterman
//! alignment score, but it will promote better translations.
//!
//! ## Shifting Eligibility
//!
//! An indel is eligible to be shifted if:
//!
//! - They must occur out of frame, starting after the first or second base of a
//!   codon.
//! - They must be a length that is a multiple of 3, so that it does not disrupt
//!   the reading frame
//!
//! These can be checked with [`CdsDeletionRange::eligible_for_shift`] and
//! [`CdsInsertionRange::eligible_for_shift`].
//!
//! For any indel eligible to be shifted, there are constraints on how far and
//! what direction it can shift. See [the deletion shifting
//! constraints](crate::outputs::Product#deletion-shifting-constraints) and [the
//! insertion shifting
//! constraints](crate::outputs::Product#insertion-shifting-constraints) for
//! more details.
//!
//! ## Shifting Direction and Amount
//!
//! An indel can be shifted to the left or to the right. When the indel appears
//! after the first base of a codon, it will either shift one base left or two
//! bases right. When the indel appears after the second base of a codon, it
//! will either shift two bases left or one base right.
//!
//! When both shift directions are possible, then codon usage statistics are
//! used to pick the most probably direction. The algorithm mirrors the logic in
//! `codonCorrectStats.pl`, albeit for coordinate math instead of strings.
//!
//! TODO: Discuss default fallback behavior for ties or non-existent stats
//!
//! ## Deletion Merging
//!
//! Two deletions may appear adjacent to each other in the product ranges,
//! caused by belonging to different exons (or caused by two deletions shifting
//! so that they are adjacent). To facilitate shifting, these adjacent deletions
//! are merged as long as they are also adjacent in reference coordinates (i.e.,
//! there is no noncoding region between them).
//!
//! ## Insertion Dropping
//!
//! Shifting may cause insertions to become at the front or end of the product,
//! or adjacent to a non-coding region between two exons. In this case, they are
//! considered part of the noncoding region and hence removed.

use crate::{
    QueryRecord,
    config::product_spec::ProductSpec,
    error,
    outputs::Product,
    ranges::{CdsDeletionRange, CdsInsertionRange, CdsMatchRange, CdsStateRange, RangeExt},
    warn,
};
use std::cmp::Ordering;

impl<'a> Product<'a> {
    /// The procedure for fixing the frames of all the eligible indels.
    ///
    /// See [the module documentation](crate::annotation::fix_frames) for more
    /// details.
    pub(crate) fn fix_frames(&mut self, query: &QueryRecord) {
        // The index of the current CdsStateRange to correct/handle
        let mut idx = 0;

        while let Some(states) = get_frame_states(idx, &mut self.product_ranges) {
            // Perform any frame fixing on range, which then returns whether to
            // advance the index or not, as well as any states to remove.
            let IdxAdjustment { advance, removal } = fix_frame(states, query, self.product_spec);

            if let Some(removal) = removal {
                // Remove the specified states, returning the resulting shift
                // that will need to be applied to idx
                let idx_shift = remove_states(self, &removal, idx);

                // Advance idx first, to avoid underflow
                if advance {
                    idx += 1;
                }

                // Apply the shift due to removed states
                idx -= idx_shift;
            } else if advance {
                // Advance idx
                idx += 1;
            }
        }
    }
}

/// The state at a given index, along with flanking states if they are in
/// bounds.
///
/// This is returned by [`get_frame_states`], and used by [`fix_frame`].
struct FrameStates<'a> {
    /// The state two to the left of the current one (`idx-2`)
    left2:   Option<&'a mut CdsStateRange>,
    /// The state left of the current one (`idx-1`)
    left1:   Option<&'a mut CdsStateRange>,
    /// The current states (`idx`)
    current: &'a mut CdsStateRange,
    /// The state right of the current one (`idx+1`)
    right1:  Option<&'a mut CdsStateRange>,
    /// The state two to the right of the current one (`idx+2`)
    right2:  Option<&'a mut CdsStateRange>,
}

/// Gets the state at `idx` within `product_ranges`, as well as the two flanking
/// states if available.
///
/// This is a helper function for [`fix_frames`]. If `None` is returned, then
/// the index is out of bounds.
///
/// [`fix_frames`]: Product::fix_frames
#[must_use]
fn get_frame_states(idx: usize, product_ranges: &mut [CdsStateRange]) -> Option<FrameStates<'_>> {
    let (left, current_and_right) = product_ranges.split_at_mut_checked(idx)?;
    let (current, right) = current_and_right.split_first_mut()?;

    let (left2, left1) = match left {
        [.., left2, left1] => (Some(left2), Some(left1)),
        [left1] => (None, Some(left1)),
        [] => (None, None),
    };

    let (right1, right2) = match right {
        [right1, right2, ..] => (Some(right1), Some(right2)),
        [right1] => (Some(right1), None),
        [] => (None, None),
    };

    Some(FrameStates {
        left2,
        left1,
        current,
        right1,
        right2,
    })
}

/// Flags indicating which of the states in [`FrameStates`] should be removed,
/// after the modifications made by [`fix_frame`].
struct StateRemoval {
    /// The state left of the current one (`idx-1`) should be removed.
    remove_left1:   bool,
    /// The current state (`idx`) should be removed.
    remove_current: bool,
    /// The state right of the current one (`idx+1`) should be removed.
    remove_right1:  bool,
    /// The state two to the right of the current one (`idx+2`) should be
    /// removed.
    remove_right2:  bool,
}

/// The return value for [`fix_frame`], indicating which states need to be
/// removed, and whether or not the index should be advanced in [`fix_frames`].
///
/// [`fix_frames`]: Product::fix_frames
struct IdxAdjustment {
    /// Whether to advance the index, or whether to rehandle the same state
    advance: bool,
    /// Any states to remove (e.g., due to becoming empty, shifting to exon
    /// boundaries, or being merged)
    removal: Option<StateRemoval>,
}

impl IdxAdjustment {
    /// Returns an [`IdxAdjustment`] that advances to the next index without
    /// removing any states.
    #[inline]
    #[must_use]
    const fn next() -> Self {
        Self {
            advance: true,
            removal: None,
        }
    }
}

/// The logic for fixing the frame of a single state within [`fix_frames`].
///
/// This will perform any mutation of the states as needed, but removal is not
/// handled by this function. Instead, [`IdxAdjustment`] is returned to indicate
/// to the caller which states need to be removed, and whether the index in
/// [`fix_frames`] should be advanced or not.
///
/// [`fix_frames`]: Product::fix_frames
#[must_use]
fn fix_frame(states: FrameStates, query: &QueryRecord, product_spec: &ProductSpec) -> IdxAdjustment {
    match states.current {
        // A match state never needs to be shifted
        CdsStateRange::M(_) => IdxAdjustment::next(),

        CdsStateRange::D(del) => match (states.left1, states.right1) {
            // A deletion flanked by match states may need to be shifted. This
            // in turn could result in empty flanking match states, so it may
            // need to be merged left if another deletion is two to the left. We
            // do not need to handle merging it to the right since our normal
            // deletion-merging branch will handle that.
            (Some(CdsStateRange::M(left_match)), Some(CdsStateRange::M(right_match))) => {
                let removal = fix_flanked_deletion(states.left2, left_match, del, right_match, query, product_spec);
                IdxAdjustment { advance: true, removal }
            }

            // A second deletion follows the current deletion, so they need to
            // be merged if possible. We then rerun the match statement on the
            // same index to allow additional merging, or to allow the deletion
            // to shift.
            (_, Some(CdsStateRange::D(right_del))) => {
                let removal = fix_adjacent_deletions(del, right_del, query, product_spec);
                IdxAdjustment {
                    // If no removal (i.e., no merging), then continue
                    advance: removal.is_none(),
                    removal,
                }
            }

            // A deletion is flanked by a deletion to the left (meaning merging
            // was not possible in the previous branches), so a warning for
            // failed shifting may be needed.
            (Some(CdsStateRange::D(_)), Some(CdsStateRange::M(_))) => {
                warn_adjacent_deletions(del, query, product_spec);
                IdxAdjustment::next()
            }

            // A deletion at the beginning of the product may require a warning
            // for failed shifting
            (None, _) => {
                warn_leading_deletion(query, del, product_spec);
                IdxAdjustment::next()
            }

            // A deletion at the end of the product may require a warning for
            // failed shifting
            (_, None) => {
                warn_trailing_deletion(query, del, product_spec);
                IdxAdjustment::next()
            }

            // A deletion is flanked by an insertion to the left or right, so a
            // warning for failed shifting may be needed.
            (Some(CdsStateRange::I(_)), _) | (_, Some(CdsStateRange::I(_))) => {
                warn_del_with_flanking_ins(del, query, product_spec);
                IdxAdjustment::next()
            }
        },

        CdsStateRange::I(ins) => match (states.left1, states.right1) {
            // An insertion flanked by match states may need to be shifted. This
            // in turn could result in empty flanking match states, so it may
            // need to be merged left or right. It could also result in the
            // indel moving to the edge of the product or exon, in which case it
            // is dropped.
            (Some(CdsStateRange::M(left_match)), Some(CdsStateRange::M(right_match))) => {
                let removal =
                    fix_flanked_insertion(states.left2, left_match, ins, right_match, states.right2, query, product_spec);
                IdxAdjustment { advance: true, removal }
            }

            // Two insertions appear adjacent to each other
            (Some(CdsStateRange::I(_)), _) | (_, Some(CdsStateRange::I(_))) => {
                warn_adjacent_insertions(query, product_spec);
                IdxAdjustment::next()
            }

            // An insertion at the beginning of the product may require a
            // warning for failed shifting
            (None, _) => {
                warn_leading_insertion(query, ins, product_spec);
                IdxAdjustment::next()
            }

            // An insertion at the end of the product may require a warning for
            // failed shifting
            (_, None) => {
                warn_trailing_insertion(query, ins, product_spec);
                IdxAdjustment::next()
            }

            // An insertion is flanked by a deletion to the left or right, so a
            // warning for failed shifting may be needed.
            (Some(CdsStateRange::D(_)), _) | (_, Some(CdsStateRange::D(_))) => {
                warn_ins_with_flanking_del(ins, query, product_spec);
                IdxAdjustment::next()
            }
        },
    }
}

/// A helper function for removing the states in a [`StateRemoval`] struct from
/// the product. The return value is the amount that must be subtracted from
/// `idx` to correct for the removed states.
#[must_use]
fn remove_states(product: &mut Product, removal: &StateRemoval, idx: usize) -> usize {
    let StateRemoval {
        remove_left1,
        remove_current,
        remove_right1,
        remove_right2,
    } = *removal;

    // Remove states from right to left
    if remove_right2 {
        product.product_ranges.remove(idx + 2);
    }
    if remove_right1 {
        product.product_ranges.remove(idx + 1);
    }
    if remove_current {
        product.product_ranges.remove(idx);
    }
    if remove_left1 {
        product.product_ranges.remove(idx - 1);
    }

    (remove_current as usize) + (remove_left1 as usize)
}

/// Performs the shifting logic for [`fix_frames`] to a deletion flanked by
/// match states.
///
/// This handles checking for shifting eligibility, performing the shift,
/// dropping empty match states, and merging the deletion left if needed. See
/// [the module documentation](crate::annotation::fix_frames) for more details.
///
/// [`fix_frames`]: Product::fix_frames
#[must_use]
fn fix_flanked_deletion(
    left2: Option<&mut CdsStateRange>, left_match: &mut CdsMatchRange, del: &mut CdsDeletionRange,
    right_match: &mut CdsMatchRange, query: &QueryRecord, product_spec: &ProductSpec,
) -> Option<StateRemoval> {
    let shift_dir = pick_deletion_shift(left_match, del, right_match, query, product_spec)?;

    apply_deletion_shift(shift_dir, left_match, del, right_match);

    // Check whether a left shift caused left_match to become empty, or a right
    // shift caused right_match to become empty (in which case it is removed)
    let (remove_left_match, remove_right_match) = if left_match.is_empty() {
        (true, false)
    } else if right_match.is_empty() {
        (false, true)
    } else {
        (false, false)
    };

    // Check whether the current deletion needs to be merged left, in which case
    // it is dropped
    let remove_del = if remove_left_match
        && let Some(CdsStateRange::D(left_del)) = left2
        && is_valid_del_merge(left_del, del, product_spec)
    {
        left_del.cds_range.end = del.cds_range.end;
        true
    } else {
        false
    };

    Some(StateRemoval {
        remove_left1:   remove_left_match,
        remove_current: remove_del,
        remove_right1:  remove_right_match,
        remove_right2:  false,
    })
}

/// Performs the shifting logic for [`fix_frames`] to an insertion flanked by
/// match states.
///
/// This handles checking for shifting eligibility, performing the shift,
/// dropping empty match states, dropping the insertion if it shifted to
/// product/exon boundaries, merging the insertion left or right if needed. See
/// [the module documentation](crate::annotation::fix_frames) for more details.
///
/// [`fix_frames`]: Product::fix_frames
fn fix_flanked_insertion(
    left2: Option<&mut CdsStateRange>, left_match: &mut CdsMatchRange, ins: &mut CdsInsertionRange,
    right_match: &mut CdsMatchRange, right2: Option<&mut CdsStateRange>, query: &QueryRecord, product_spec: &ProductSpec,
) -> Option<StateRemoval> {
    let shift_dir = pick_insertion_shift(left_match, ins, right_match, query, product_spec)?;

    apply_insertion_shift(shift_dir, left_match, ins, right_match);

    // Flags which will be set if particular states need to be removed
    let mut remove_ins = false;
    let mut remove_left_match = false;
    let mut remove_right_match = false;
    let mut remove_right2 = false;

    // Check whether the current insertion is now adjacent to a non-coding
    // region
    let ins_adj_noncoding = product_spec
        .exons
        .noncoding_regions
        .iter()
        .any(|noncoding| noncoding.cds_index == ins.cds_index);

    if ins_adj_noncoding {
        remove_ins = true;
    }

    // Check whether a left shift caused left_match to become empty, or a right
    // shift caused right_match to become empty
    if left_match.is_empty() {
        // Remove the empty match state
        remove_left_match = true;

        // Handle the consequences of removing the empty match state
        match left2 {
            // Merge the now-adjacent insertions
            Some(CdsStateRange::I(left_ins)) if !remove_ins => {
                validate_ins_merge(left_ins, ins);
                left_ins.query_range.end = ins.query_range.end;
                remove_ins = true;
            }

            // Remove insertion if it shifted to the beginning of the product
            None => {
                remove_ins = true;
            }

            Some(_) => {}
        }
    } else if right_match.is_empty() {
        // Remove the empty match state
        remove_right_match = true;

        // Handle the consequences of removing the empty match state
        match right2 {
            // Merge the now-adjacent insertions
            Some(CdsStateRange::I(right_ins)) if !remove_ins => {
                validate_ins_merge(ins, right_ins);
                ins.query_range.end = right_ins.query_range.end;
                remove_right2 = true;
            }

            // Remove insertion if it shifted to the end of the product
            None => remove_ins = true,

            Some(_) => {}
        }
    }

    Some(StateRemoval {
        remove_left1: remove_left_match,
        remove_current: remove_ins,
        remove_right1: remove_right_match,
        remove_right2,
    })
}

/// Performs the deletion merging logic for [`fix_frames`] to two adjacent
/// deletions, specifically a deletion occurring to the right of the current
/// deletion.
///
/// This handles merging the deletions if eligible and issuing a warning with
/// [`warn_adjacent_deletions`] if not eligible.
///
/// [`fix_frames`]: Product::fix_frames
fn fix_adjacent_deletions(
    del: &mut CdsDeletionRange, right_del: &mut CdsDeletionRange, query: &QueryRecord, product_spec: &ProductSpec,
) -> Option<StateRemoval> {
    // Validity: product_ranges partitions the aligned-against
    // coding sequence, so adjacent deletions in product_ranges are
    // also adjacent in the coding sequence
    debug_assert_eq!(del.cds_range.end, right_del.cds_range.start);

    // Validity: TODO
    if is_valid_del_merge(del, right_del, product_spec) {
        // Extend range to encompass right_range
        del.cds_range.end = right_del.cds_range.end;

        return Some(StateRemoval {
            remove_left1:   false,
            remove_current: false,
            remove_right1:  true,
            remove_right2:  false,
        });
    }

    warn_adjacent_deletions(del, query, product_spec);
    None
}

/// Warns about failure to shift a deletion if another adjacent deletion is
/// present and the current deletion is eligible for shifting.
///
/// Notably, if the current deletion has another deletion to its left, this
/// means the two deletions are not able to be merged (since otherwise they
/// would've been merged when `idx` was `idx-1` in [`fix_frames`]).
///
/// [`fix_frames`]: Product::fix_frames
fn warn_adjacent_deletions(del: &CdsDeletionRange, query: &QueryRecord, product_spec: &ProductSpec) {
    // Issue warning if applicable that shift cannot be
    // performed
    if del.eligible_for_shift() {
        warn!(
            "Two deletions within adjacent exons occurred. Because they were separated by a non-coding region, frame correction cannot be applied. Consider manually adjusting the output alignment.\nQuery: {query_id}\nProduct: {product}",
            query_id = query.id,
            product = product_spec.name
        )
    }
}

/// Warns about two adjacent insertions, which should never occur in the input
/// product ranges.
fn warn_adjacent_insertions(query: &QueryRecord, product_spec: &ProductSpec) {
    error!(
        "Two insertions occurred adjacent to each other.\nQuery: {query_id}\nProduct: {product}",
        query_id = query.id,
        product = product_spec.name
    );
}

/// Warns about a deletion with a flanking insertion, if the deletion is
/// eligible to be shifted.
fn warn_del_with_flanking_ins(del: &CdsDeletionRange, query: &QueryRecord, product_spec: &ProductSpec) {
    if del.eligible_for_shift() {
        warn!(
            "A deletion occurred adjacent to an insertion, so frame correction cannot be applied. Consider manually adjusting the output alignment.\nQuery: {query_id}\nProduct: {product}",
            query_id = query.id,
            product = product_spec.name
        );
    }
}

/// Warns about an insertion with a flanking deletion, if the insertion is
/// eligible to be shifted.
fn warn_ins_with_flanking_del(ins: &CdsInsertionRange, query: &QueryRecord, product_spec: &ProductSpec) {
    if ins.eligible_for_shift() {
        warn!(
            "An insertion occurred adjacent to a deletion, so frame correction cannot be applied. Consider manually adjusting the output alignment.\nQuery: {query_id}\nProduct: {product}",
            query_id = query.id,
            product = product_spec.name
        );
    }
}

/// Warns about an insertion at the beginning of a product, if the insertion is
/// eligible to be shifted.
fn warn_leading_insertion(query: &QueryRecord, ins: &CdsInsertionRange, product_spec: &ProductSpec) {
    if ins.eligible_for_shift() {
        warn!(
            "An insertion appeared at the start of the product alignment, so frame correction cannot be applied. Consider manually adjusting the output alignment.\nQuery: {query_id}\nProduct: {product}",
            query_id = query.id,
            product = product_spec.name
        );
    }
}

/// Warns about a deletion at the beginning of a product, if the deletion is
/// eligible to be shifted.
fn warn_leading_deletion(query: &QueryRecord, del: &CdsDeletionRange, product_spec: &ProductSpec) {
    if del.eligible_for_shift() {
        warn!(
            "A deletion appeared at the start of the product alignment, so frame correction cannot be applied. Consider manually adjusting the output alignment.\nQuery: {query_id}\nProduct: {product}",
            query_id = query.id,
            product = product_spec.name
        );
    }
}

/// Warns about an insertion at the end of a product, if the insertion is
/// eligible to be shifted.
fn warn_trailing_insertion(query: &QueryRecord, ins: &CdsInsertionRange, product_spec: &ProductSpec) {
    if ins.eligible_for_shift() {
        warn!(
            "An insertion appeared at the end of the product alignment, so frame correction cannot be applied. Consider manually adjusting the output alignment.\nQuery: {query_id}\nProduct: {product}",
            query_id = query.id,
            product = product_spec.name
        );
    }
}

/// Warns about a deletion at the end of a product, if the deletion is eligible
/// to be shifted.
fn warn_trailing_deletion(query: &QueryRecord, del: &CdsDeletionRange, product_spec: &ProductSpec) {
    if del.eligible_for_shift() {
        warn!(
            "A deletion appeared at the end of the product alignment, so frame correction cannot be applied. Consider manually adjusting the output alignment.\nQuery: {query_id}\nProduct: {product}",
            query_id = query.id,
            product = product_spec.name
        );
    }
}

/// Determines whether two adjacent deletions are eligible to be merged.
///
/// To be merged, they must appear on the same side of every noncoding region,
/// so that the merged deletion does not span such a region. Were a deletion to
/// span a noncoding region, then shifting it could cause bases in the noncoding
/// region to be moved into flanking match states.
///
/// ## Validity
///
/// The two deletions passed must be adjacent in CDS coordinates. They also
/// should not span any noncoding regions.
fn is_valid_del_merge(del1: &CdsDeletionRange, del2: &CdsDeletionRange, product_spec: &ProductSpec) -> bool {
    // Validity: product_ranges partitions the aligned-against
    // coding sequence, so adjacent deletions in product_ranges are
    // also adjacent in the coding sequence
    debug_assert_eq!(del1.cds_range.end, del2.cds_range.start);

    // For each noncoding region, the deletions must both come before it, or
    // must both come after it. Neither deletion should contain a noncoding
    // region
    product_spec
        .exons
        .noncoding_regions
        .iter()
        .all(|noncoding| del1.cds_range.cmp_ins(&noncoding.cds_index) == del2.cds_range.cmp_ins(&noncoding.cds_index))
}

/// Performs checking of two adjacent insertions to confirm they are eligible to
/// be merged, issuing errors if not.
fn validate_ins_merge(ins1: &CdsInsertionRange, ins2: &CdsInsertionRange) {
    if ins1.cds_index != ins2.cds_index {
        error!("Two adjacent insertions had different CDS indices");
    } else if ins1.query_range.end != ins2.query_range.start {
        // TODO: We need to make sure to drop insertions shifting to edges
        // of exons or products
        error!("Two adjacent insertions were not adjacent in the query sequence");
    }
}

/// The chosen direction to shift an out-of-frame indel.
enum ShiftDir {
    /// The indel should be shifted left (causing 1-2 match states to move
    /// right).
    Left,
    /// The indel should be shifted right (causing 1-2 match states to move
    /// left).
    Right,
}

/// Picks the direction to shift an insertion if it is out of frame, or
/// returns `None` if the no shift is needed or allowed.
///
/// Formally, this method looks at the constraints for how the insertion is
/// allowed to shift, and evaluates whether it is eligible to shift left or
/// right. If both are allowed, then codon usage statistics are used to
/// determine the direction. If only one is allowed, then that shift
/// direction is used. If neither are allowed, `None` is returned.
///
/// ## Insertion Shifting Constraints
///
/// Currently, if an insertion is adjacent to a non-coding region, then no
/// shift is performed. This could be modified TODO.
///
/// The magnitude of the shift is limited by the number of flanking match
/// states in the direction of the shift, since bases will need to be moved
/// from the match states into the insert state.
///
/// If the insertion is strictly contained inside a region where adjacent
/// exons overlap, then the insertion may not shift (since it actually
/// appears twice in the CDS).
///
/// Similarly, an insertion cannot shift into one of the overlapping
/// regions.
fn pick_insertion_shift(
    left_match: &CdsMatchRange, ins: &CdsInsertionRange, right_match: &CdsMatchRange, query: &QueryRecord,
    product_spec: &ProductSpec,
) -> Option<ShiftDir> {
    let frameshift = ins.frameshift();

    // Check shift eligibility, equivalent to eligible_for_shift but
    // allowing us to retain the frameshift.
    if frameshift == 0 || !ins.len().is_multiple_of(3) {
        return None;
    }

    // TODO: Reevaluate whether forced shift is possible. No decision was
    // ever made for this, but I am writing this to mirror deletion case
    //
    // Check adjacency to non-coding region, which prevents shift
    if left_match.cds_range.end != right_match.cds_range.start {
        warn!(
            "An insertion adjacent to a non-coding region was found and will not be shifted.\nQuery: {query_id}\nProduct: {product}",
            query_id = query.id,
            product = product_spec.name
        );
        return None;
    }

    // Get the maximum left shift or right shift allowed, based on
    // insufficient flanking state, nearby overlapping exons, or spanned
    // overlapping exons (or 2 if no constraint).
    let (max_left_shift, max_right_shift) = {
        // A shift will never be by more than 2 in all cases
        let mut max_left_shift = 2;
        let mut max_right_shift = 2;

        // Limit max_left_shift by the number of flanking match states left
        // of it
        max_left_shift = max_left_shift.min(left_match.len());

        // Limit max_right_shift by the number of flanking match states
        // right of it
        max_right_shift = max_right_shift.min(right_match.len());

        // Check if the insertion is inside an overlap region
        let inside_overlap_region = product_spec.exons.overlapped_regions.iter().any(|overlap| {
            // TODO: Same logic as in InsertionRange::intersect_exon, so
            // abstract
            overlap.cds_range().start < ins.cds_index.right() && ins.cds_index.right() < overlap.cds_range().end
        });

        // If we are inside an overlap region, no shift can occur since the
        // insertion is duplicated in the CDS
        if inside_overlap_region {
            return None;
        }

        // Get the number of bases between the insertion and the closest
        // overlapping region left of the insertion. Use next_back to
        // priorize the rightmost overlapping region. If an overlap region
        // is directly left of the insertion, then
        // ins.cds_index.right()-overlap.cds_range().end == 0, hence the use
        // of right() instead of left().
        let left_adj_overlap_dist = product_spec
            .exons
            .overlapped_regions
            .iter()
            .filter_map(|overlap| ins.cds_index.right().checked_sub(overlap.cds_range().end))
            .next_back();

        // Limit max_left_shift based on the closest overlapping region left
        // of the insertion.
        if let Some(left_adj_overlap_dist) = left_adj_overlap_dist {
            max_left_shift = max_left_shift.min(left_adj_overlap_dist);
        }

        // Get the number of bases between the insertion and the closest
        // overlapping region right of the insertion. Use next to prioritize
        // the leftmost overlapping region. If an overlap region is directly
        // right of the insertion, then
        // overlap.cds_range().start-ins.cds_index.right() == 0, hence the
        // use of right() instead of left()
        let right_adj_overlap_dist = product_spec
            .exons
            .overlapped_regions
            .iter()
            .filter_map(|overlap| overlap.cds_range().start.checked_sub(ins.cds_index.right()))
            .next();

        // Limit max_right_shift based on the closest overlapping region
        // right of the insertion.
        if let Some(right_adj_overlap_dist) = right_adj_overlap_dist {
            max_right_shift = max_right_shift.min(right_adj_overlap_dist);
        }

        (max_left_shift, max_right_shift)
    };

    // Check whether shifts are possible given the frame
    let (can_shift_left, can_shift_right) = if frameshift == 1 {
        (max_left_shift >= 1, max_right_shift >= 2)
    } else if frameshift == 2 {
        (max_left_shift >= 2, max_right_shift >= 1)
    } else {
        return None;
    };

    // Check for forced or impossible shifts
    match (can_shift_left, can_shift_right) {
        (false, false) => {
            warn!(
                "An insertion could not be shifted due to insufficient flanking match state or the presence of overlapping exons.\nQuery: {query_id}\nProduct: {product}",
                query_id = query.id,
                product = product_spec.name
            );
            return None;
        }
        (true, false) => {
            warn!(
                "An insertion is forcibly being shifted left (e.g., due to insufficient flanking match state or the presence of overlapping exons).\nQuery: {query_id}\nProduct: {product}",
                query_id = query.id,
                product = product_spec.name
            );
            return Some(ShiftDir::Left);
        }
        (false, true) => {
            warn!(
                "An insertion is forcibly being shifted right (e.g., due to insufficient flanking match state or the presence of overlapping exons).\nQuery: {query_id}\nProduct: {product}",
                query_id = query.id,
                product = product_spec.name
            );
            return Some(ShiftDir::Right);
        }
        (true, true) => {}
    }

    let codon_position = ins.cds_index.to_aa_idx().right_pos();
    let insert_seq = &query.nucleotides[ins.query_range.clone()];

    if frameshift == 2 {
        let [i1, .., i2, i3] = *insert_seq else {
            return None;
        };

        let cp1 = query.nucleotides[left_match.query_range.end - 2];
        let cp2 = query.nucleotides[left_match.query_range.end - 1];
        let cp3 = query.nucleotides[right_match.query_range.start];

        // Codon formed by shifting insertion right 1
        let a2r1 = build_discontiguous_codon(cp1, cp2, i1);
        // Codon formed by shifting insertion left 2
        let a2l2 = build_discontiguous_codon(i2, i3, cp3);

        // TODO: What about N vs X
        // Validity: The codons are uppercase because they are derived from
        // bases in query
        if product_spec.codon_left_ge_right(a2r1, a2l2, codon_position as u32) {
            Some(ShiftDir::Right)
        } else {
            Some(ShiftDir::Left)
        }
    } else if frameshift == 1 {
        let [i1, i2, .., i3] = *insert_seq else {
            return None;
        };

        let cp1 = query.nucleotides[left_match.query_range.end - 1];
        let cp2 = query.nucleotides[right_match.query_range.start];
        let cp3 = query.nucleotides[right_match.query_range.start + 1];

        // Codon formed by shifting insertion right 2
        let a1r2 = build_discontiguous_codon(cp1, i1, i2);
        // Codon formed by shifting insertion left 1
        let a1l1 = build_discontiguous_codon(i3, cp2, cp3);

        // Validity: build_discontiguous_codon ensures validity requirements
        // are met
        if product_spec.codon_left_ge_right(a1l1, a1r2, codon_position as u32) {
            Some(ShiftDir::Left)
        } else {
            Some(ShiftDir::Right)
        }
    } else {
        None
    }
}

/// Picks the direction to shift a deletion if it is out of frame, or
/// returns `None` if the no shift is needed or allowed.
///
/// Formally, this method looks at the constraints for how the deletion is
/// allowed to shift, and evaluates whether it is eligible to shift left or
/// right. If both are allowed, then codon usage statistics are used to
/// determine the direction. If only one is allowed, then that shift
/// direction is used. If neither are allowed, `None` is returned.
///
/// ## Deletion Shifting Constraints
///
/// Currently, if a deletion is adjacent to a non-coding region, then no
/// shift is performed. This could be modified TODO.
///
/// The magnitude of the shift is limited by the number of flanking match
/// states in the direction of the shift, since bases will need to be moved
/// from the match states into the delete state.
///
/// If the deletion crosses over a region where adjacent exons overlap, then
/// the deletion may not shift in a way that would cause the deletion to no
/// longer fully span that region.
///
/// Similarly, a deletion cannot shift in a way that causes it to intersect
/// one of the overlapping regions when it did not before.
fn pick_deletion_shift(
    left_match: &CdsMatchRange, del: &CdsDeletionRange, right_match: &CdsMatchRange, query: &QueryRecord,
    product_spec: &ProductSpec,
) -> Option<ShiftDir> {
    let frameshift = del.frameshift();

    // Check shift eligibility, equivalent to eligible_for_shift but
    // allowing us to retain the frameshift.
    if frameshift == 0 || !del.len().is_multiple_of(3) {
        return None;
    }

    // TODO: Reevaluate whether forced shift is possible for adjacent case.
    //
    // Check adjacency to non-coding region, which prevents shift
    if left_match.query_range.end != right_match.query_range.start {
        warn!(
            "A deletion adjacent to or crossing into a non-coding region was found and will not be shifted.\nQuery: {query_id}\nProduct: {product}",
            query_id = query.id,
            product = product_spec.name
        );
        return None;
    }

    // Get the maximum left shift or right shift allowed, based on
    // insufficient flanking state, nearby overlapping exons, or spanned
    // overlapping exons (or 2 if no constraint).
    let (max_left_shift, max_right_shift) = {
        // A shift will never be by more than 2 in all cases
        let mut max_left_shift = 2;
        let mut max_right_shift = 2;

        // Limit max_left_shift by the number of flanking match states left
        // of it
        max_left_shift = max_left_shift.min(left_match.len());

        // Limit max_right_shift by the number of flanking match states
        // right of it
        max_right_shift = max_right_shift.min(right_match.len());

        // Limit max_left_shift and max_right_shift due to the deletion
        // spanning overlapping exons
        let mut spanned_overlaps =
            product_spec.exons.overlapped_regions.iter().filter(|overlap| {
                del.cds_range.overlaps(&overlap.cds_range1) || del.cds_range.overlaps(&overlap.cds_range2)
            });
        if let Some(first) = spanned_overlaps.next() {
            let last = spanned_overlaps.next_back().unwrap_or(first);

            let overlap_cds_start = first.cds_range1.start;
            let overlap_cds_end = last.cds_range2.end;

            let deletion_before_overlap = overlap_cds_start - del.cds_range.start;
            let deletion_after_overlap = del.cds_range.end - overlap_cds_end;

            max_left_shift = max_left_shift.min(deletion_before_overlap);
            max_right_shift = max_right_shift.min(deletion_after_overlap);
        }

        // Get the number of bases between the deletion and the closest
        // overlapping region left of the deletion. Use next_back to
        // priorize the rightmost overlapping region.
        let left_adj_overlap_dist = product_spec
            .exons
            .overlapped_regions
            .iter()
            .filter_map(|overlap| del.cds_range.start.checked_sub(overlap.cds_range().end))
            .next_back();

        // Limit max_left_shift based on the closest overlapping region left
        // of the deletion.
        if let Some(left_adj_overlap_dist) = left_adj_overlap_dist {
            max_left_shift = max_left_shift.min(left_adj_overlap_dist);
        }

        // Get the number of bases between the deletion and the closest
        // overlapping region right of the deletion. Use next to prioritize
        // the leftmost overlapping region.
        let right_adj_overlap_dist = product_spec
            .exons
            .overlapped_regions
            .iter()
            .filter_map(|overlap| overlap.cds_range().start.checked_sub(del.cds_range.end))
            .next();

        // Limit max_right_shift based on the closest overlapping region
        // right of the deletion.
        if let Some(right_adj_overlap_dist) = right_adj_overlap_dist {
            max_right_shift = max_right_shift.min(right_adj_overlap_dist);
        }

        (max_left_shift, max_right_shift)
    };

    // Check whether shifts are possible given the frame
    let (can_shift_left, can_shift_right) = if frameshift == 1 {
        (max_left_shift >= 1, max_right_shift >= 2)
    } else if frameshift == 2 {
        (max_left_shift >= 2, max_right_shift >= 1)
    } else {
        return None;
    };

    // Check for forced or impossible shifts
    match (can_shift_left, can_shift_right) {
        (false, false) => {
            warn!(
                "A deletion could not be shifted due to insufficient flanking match state or the presence of overlapping exons.\nQuery: {query_id}\nProduct: {product}",
                query_id = query.id,
                product = product_spec.name
            );
            return None;
        }
        (true, false) => {
            warn!(
                "A deletion is forcibly being shifted left (e.g., due to insufficient flanking match state or the presence of overlapping exons).\nQuery: {query_id}\nProduct: {product}",
                query_id = query.id,
                product = product_spec.name
            );
            return Some(ShiftDir::Left);
        }
        (false, true) => {
            warn!(
                "A deletion is forcibly being shifted right (e.g., due to insufficient flanking match state or the presence of overlapping exons).\nQuery: {query_id}\nProduct: {product}",
                query_id = query.id,
                product = product_spec.name
            );
            return Some(ShiftDir::Right);
        }
        (true, true) => {}
    }

    // Both shifts are possible, so get the codon positions (1-based) at the
    // boundaries of the deletion for table checking
    let pos_left = (del.cds_range.start / 3) + 1;
    let pos_right = ((del.cds_range.end - 1) / 3) + 1;

    if frameshift == 1 {
        // We need to shift the deletion left by 1 or right by 2

        let cp1 = query.nucleotides[left_match.query_range.end - 1];
        let cp2 = query.nucleotides[right_match.query_range.start];
        let cp3 = query.nucleotides[right_match.query_range.start + 1];

        let pivot = build_discontiguous_codon(cp1, cp2, cp3);

        // TODO: Likely want to compare with Ordering::is_gt instead, so
        // that ties resolve as right shift.

        // By default, we shift the deletion left for frame 1 (causing 1
        // match to shift right, rather than 2). Only if there is evidence
        // for shifting the deletion to right do we do that. pos_left being
        // better means shifting matches left is better, which means
        // shifting deletion right is better. Validity: The codon is
        // uppercase because it is derived from bases in query
        let shift_del_right = product_spec
            .compare_codon_positions(pos_left as u32, pos_right as u32, pivot)
            .is_some_and(Ordering::is_ge);

        if shift_del_right {
            Some(ShiftDir::Right)
        } else {
            Some(ShiftDir::Left)
            // Shift the deletion left by 1, causing 1 match to move right
        }
    } else if frameshift == 2 {
        let cp1 = query.nucleotides[left_match.query_range.end - 2];
        let cp2 = query.nucleotides[left_match.query_range.end - 1];
        let cp3 = query.nucleotides[right_match.query_range.start];

        let pivot = build_discontiguous_codon(cp1, cp2, cp3);

        // By default, we shift the deletion right for frame 2 (causing 1
        // match to shift left, rather than 2). Only if there is evidence
        // for a left shift do we do that. pos_right being better means
        // shifting matches right is better, which means shifting deletion
        // left is better. Validity: The codon is uppercase because it is
        // derived from bases in query
        let shift_del_left = product_spec
            .compare_codon_positions(pos_left as u32, pos_right as u32, pivot)
            .is_some_and(Ordering::is_lt);

        if shift_del_left {
            Some(ShiftDir::Left)
        } else {
            Some(ShiftDir::Right)
        }
    } else {
        None
    }
}

/// Given a shifting direction chosen by [`pick_insertion_shift`], applies the
/// shift to the adjacent states in the product ranges.
fn apply_insertion_shift(
    shift_dir: ShiftDir, left_match: &mut CdsMatchRange, ins: &mut CdsInsertionRange, right_match: &mut CdsMatchRange,
) {
    let frameshift = ins.frameshift();

    if frameshift == 1 {
        match shift_dir {
            ShiftDir::Left => {
                left_match.cut_end(1);
                ins.shift_left(1);
                right_match.extend_start(1);
            }
            ShiftDir::Right => {
                left_match.extend_end(2);
                ins.shift_right(2);
                right_match.cut_start(2);
            }
        }
    } else if frameshift == 2 {
        match shift_dir {
            ShiftDir::Left => {
                left_match.cut_end(2);
                ins.shift_left(2);
                right_match.extend_start(2);
            }
            ShiftDir::Right => {
                left_match.extend_end(1);
                ins.shift_right(1);
                right_match.cut_start(1);
            }
        }
    }
}

/// Given a shifting direction chosen by [`pick_deletion_shift`], applies the
/// shift to the adjacent states in the product ranges.
fn apply_deletion_shift(
    shift_dir: ShiftDir, left_match: &mut CdsMatchRange, del: &mut CdsDeletionRange, right_match: &mut CdsMatchRange,
) {
    let frameshift = del.frameshift();

    if frameshift == 1 {
        match shift_dir {
            ShiftDir::Left => {
                left_match.cut_end(1);
                del.shift_left(1);
                right_match.extend_start(1);
            }
            ShiftDir::Right => {
                left_match.extend_end(2);
                del.shift_right(2);
                right_match.cut_start(2);
            }
        }
    } else if frameshift == 2 {
        match shift_dir {
            ShiftDir::Left => {
                left_match.cut_end(2);
                del.shift_left(2);
                right_match.extend_start(2);
            }
            ShiftDir::Right => {
                left_match.extend_end(1);
                del.shift_right(1);
                right_match.cut_start(1);
            }
        }
    }
}

/// Combines three discontinuous bases into a codon, automatically converting
/// `U` to `T`.
///
/// No other bytes are altered, and case is not changed.
fn build_discontiguous_codon(b1: u8, b2: u8, b3: u8) -> [u8; 3] {
    [b1, b2, b3].map(|base| if base == b'U' { b'T' } else { base })
}

impl CdsDeletionRange {
    /// Shifts the deletion in the coding sequence to the left by `amount`.
    ///
    /// This subtracts `amount` from the start and end of the range.
    fn shift_left(&mut self, amount: usize) {
        self.cds_range = self.cds_range.sub(amount);
    }

    /// Shifts the deletion in the coding sequence to the right by `amount`.
    ///
    /// This adds `amount` to the start and end of the range.
    fn shift_right(&mut self, amount: usize) {
        self.cds_range = self.cds_range.add(amount);
    }

    /// The frame shift of the deletion.
    fn frameshift(&self) -> usize {
        self.cds_range.start % 3
    }

    /// Returns whether a deletion is eligible to be shifted based solely on its
    /// frame and length.
    ///
    /// The length of the deletion in CDS coordinates must be a multiple of
    /// three, and frame must be non-zero.
    ///
    /// Note that the result of this can vary before and after a deletion is
    /// merged with adjacent deletions.
    fn eligible_for_shift(&self) -> bool {
        self.frameshift() != 0 && self.len().is_multiple_of(3)
    }
}

impl CdsInsertionRange {
    fn shift_left(&mut self, amount: usize) {
        *self.cds_index.right_mut() -= amount;
        self.query_range = self.query_range.start - amount..self.query_range.end - amount;
    }

    fn shift_right(&mut self, amount: usize) {
        *self.cds_index.right_mut() += amount;
        self.query_range = self.query_range.start + amount..self.query_range.end + amount;
    }

    /// The frame shift of the insertion.
    fn frameshift(&self) -> usize {
        self.cds_index.codon_shift()
    }

    /// Returns whether an insertion is eligible to be shifted based solely on
    /// its frame and length.
    ///
    /// The length of the insertion in CDS coordinates must be a multiple of
    /// three, and frame must be non-zero.
    fn eligible_for_shift(&self) -> bool {
        self.frameshift() != 0 && self.len().is_multiple_of(3)
    }
}
