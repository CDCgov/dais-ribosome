use crate::{
    AnnotationModule,
    annotation::rewrites::get_states::{StateWithFlanking, get_state_with_flanking},
    config::{annotation_module::ReferenceGroup, rewrite::RewriteRanges},
    ranges::{MatchRange, StateRange},
};
use zoe::alignment::AlignmentStates;

impl<'a> AnnotationModule<'a> {
    /// Rewrite deletions in the alignment per the experimental rewrite rules.
    ///
    /// Rewriting will forcibly place a deletion to be in a new location,
    /// regardless of score. Rewriting only occurs if:
    ///
    /// - The `from` range is a deletion with exactly the specified length (it
    ///   cannot be longer than the `from` range)
    /// - The `to` range must be included in the alignment
    /// - Positions in the `to` range that are not also in `from` must be
    ///   matches
    /// - Any positions between the two ranges must also be matches
    /// - At least one additional match state must be present adjacent to the
    ///   `to` range to prevent pathological cases
    ///
    /// Multiple rules are applied from left to right based on the `from` range.
    // TODO: Fix this to be more efficient when no deletions are present. This
    // could potentially be to change the order in which rules are applied so
    // that only a single pass is needed
    pub(crate) fn rule_rewrite_dels(
        &self, genome_aln_states: &mut Vec<StateRange>, alignment: &mut AlignmentStates, reference: &ReferenceGroup,
        query_len: usize,
    ) {
        // Track whether any rewriting rule were applied, in which case we
        // recompute the Zoe alignment states
        let mut any_rewritten = false;

        // Apply each rule in the order
        for rewrite_del in &reference.rewrite_dels {
            let rewritten = rewrite_del.rewrite_deletion(genome_aln_states);
            any_rewritten |= rewritten;
        }

        if any_rewritten {
            *alignment = StateRange::state_ranges_to_alignment_states(genome_aln_states, query_len);
        }
    }
}

/// The direction that a deletion reposition rule will shift the deletion.
enum RepositionDir {
    /// The `to` range is left of the `from` range.
    Left,
    /// The `to` range is right of the `from` range.
    Right,
}

impl RewriteRanges {
    /// Returns the direction of the rewrite rule.
    fn direction(&self) -> RepositionDir {
        // Validity: the ranges are the same length but not equal, so comparing
        // start points is sufficient (and no equality case is necessary)
        if self.to.start > self.from.start {
            RepositionDir::Right
        } else {
            RepositionDir::Left
        }
    }

    /// Returns the distance (number of nucleotides) that a deletion will shift
    /// under the given rule.
    fn distance(&self) -> usize {
        // Validity: to ≠ from AND to.len = from.len, so subtracting the
        // starting indices is sufficient
        self.to.start.abs_diff(self.from.start)
    }

    /// Applies the rewrite rule to the given alignment for deletions, returning
    /// `true` if the rule was applied (and the alignment mutated) or `false` if
    /// it was not applicable.
    fn rewrite_deletion(&self, genome_aln_states: &mut Vec<StateRange>) -> bool {
        // Find a deletion whose range equals `from`
        let Some(del_idx) = genome_aln_states
            .iter()
            .position(|state| matches!(state, StateRange::D(del) if del.ref_range == self.from))
        else {
            return false;
        };

        // Get the flanking states
        let Some(StateWithFlanking {
            left2: _,
            left1,
            current: StateRange::D(del),
            right1,
            right2: _,
        }) = get_state_with_flanking(del_idx, genome_aln_states)
        else {
            // This should be unreachable, since del_idx is known to be in
            // bounds and correspond to a deletion
            return false;
        };

        match self.direction() {
            RepositionDir::Left => {
                // Left of the deletion must be a match state
                let Some(StateRange::M(left1)) = left1 else { return false };

                // The match state must have sufficient length to cover the new
                // `to` range
                let distance = self.distance();
                // We require a strictly greater match length to prevent
                // pathological cases such as shifting to edge of alignment,
                // shifting adjacent to another deletion, etc. TODO: Convert to
                // strict inequality.
                if left1.len() <= distance {
                    return false;
                }

                // Shrink the match state and shift the deletion
                left1.cut_end(distance);
                del.shift_left(distance);

                if let Some(StateRange::M(right1)) = right1 {
                    // Extend the right match state
                    right1.extend_start(distance);
                } else {
                    // Validity: this will be in bounds for the query and
                    // reference based on the cut_end and shift_left calls above
                    // (this is effectively undoing that on one side)
                    let query_range = left1.query_range.end..left1.query_range.end + distance;
                    let ref_range = del.ref_range.end..del.ref_range.end + distance;

                    // Insert after the deletion
                    genome_aln_states.insert(del_idx + 1, StateRange::M(MatchRange { query_range, ref_range }));
                }
            }
            RepositionDir::Right => {
                // Right of the deletion must be a match state
                let Some(StateRange::M(right1)) = right1 else { return false };

                // The match state must have sufficient length to cover the new
                // `to` range
                let distance = self.distance();
                // We require a strictly greater match length to prevent
                // pathological cases such as shifting to edge of alignment,
                // shifting adjacent to another deletion, etc. TODO: Convert to
                // strict inequality.
                if right1.len() <= distance {
                    return false;
                }

                // Shrink the match state and shift the deletion
                right1.cut_start(distance);
                del.shift_right(distance);

                if let Some(StateRange::M(left1)) = left1 {
                    // Extend the left match state
                    left1.extend_end(distance);
                } else {
                    // Validity: this will be in bounds for the query and
                    // reference based on the cut_start and shift_right calls
                    // above (this is effectively undoing that on one side)
                    let query_range = right1.query_range.start - distance..right1.query_range.start;
                    let ref_range = del.ref_range.start - distance..del.ref_range.start;

                    // Insert at the index of the deletion, so that it appears
                    // before
                    genome_aln_states.insert(del_idx, StateRange::M(MatchRange { query_range, ref_range }));
                }
            }
        }

        true
    }
}
