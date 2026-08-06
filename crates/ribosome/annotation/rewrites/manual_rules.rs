use std::cmp::Ordering;

use crate::{
    AnnotationModule,
    annotation::rewrites::get_states::{IdxAdjustment, StateVecEdits, StateWithFlanking, rewrite},
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

        let mut rewrite_del_rule = reference.rewrite_dels.as_slice();

        rewrite(genome_aln_states, |states| {
            apply_rewrite_dels_at_state(states, &mut rewrite_del_rule, &mut any_rewritten)
        });

        if any_rewritten {
            *alignment = StateRange::state_ranges_to_alignment_states(genome_aln_states, query_len);
        }
    }
}

fn apply_rewrite_dels_at_state(
    states: StateWithFlanking<StateRange>, rewrite_del_rule: &mut &[RewriteRanges], any_rewritten: &mut bool,
) -> IdxAdjustment<StateRange> {
    // All branches of this function currently return IdxAdjustment with advance
    // as true. This means that the current deletion will be advanced past, so
    // rule chaining (applying two rules to the same deletion) does not occur

    // Ensure the current state is a deletion, or skip to next state
    let StateWithFlanking {
        left2: _,
        left1,
        current: StateRange::D(del),
        right1,
        right2: _,
    } = states
    else {
        return IdxAdjustment::next();
    };

    // Loop until we find the next applicable rule
    let rule = loop {
        let Some((rule, rest)) = rewrite_del_rule.split_first() else {
            return IdxAdjustment::next();
        };

        let ordering = rule
            .from
            .start
            .cmp(&del.ref_range.start)
            .then_with(|| rule.from.end.cmp(&del.ref_range.end));

        match ordering {
            // The rule occurs for an earlier range, so remove the rule
            // and continue loop
            Ordering::Less => *rewrite_del_rule = rest,

            // An applicable rule was found, so return the rule
            Ordering::Equal => {
                *rewrite_del_rule = rest;
                break rule;
            }

            // No rules apply, so skip to next state
            Ordering::Greater => return IdxAdjustment::next(),
        }
    };

    match rule.direction() {
        RepositionDir::Left => {
            // Left of the deletion must be a match state
            let Some(StateRange::M(left)) = left1 else {
                return IdxAdjustment::next();
            };

            // The match state must have sufficient length to cover the new
            // `to` range
            let distance = rule.distance();
            // We require a strictly greater match length to prevent
            // pathological cases such as shifting to edge of alignment,
            // shifting adjacent to another deletion, etc. TODO: Convert to
            // strict inequality.
            if left.len() <= distance {
                return IdxAdjustment::next();
            }

            *any_rewritten = true;

            // Shrink the match state and shift the deletion
            left.cut_end(distance);
            del.shift_left(distance);

            if let Some(StateRange::M(right)) = right1 {
                // Extend the right match state
                right.extend_start(distance);

                IdxAdjustment::next()
            } else {
                // Validity: this will be in bounds for the query and
                // reference based on the cut_end and shift_left calls above
                // (this is effectively undoing that on one side)
                let query_range = left.query_range.end..left.query_range.end + distance;
                let ref_range = del.ref_range.end..del.ref_range.end + distance;

                let new_match_state = StateRange::M(MatchRange { query_range, ref_range });

                // Insert after the deletion
                IdxAdjustment {
                    advance: true,
                    edits:   Some(StateVecEdits {
                        insert_right: Some(new_match_state),
                        ..Default::default()
                    }),
                }
            }
        }
        RepositionDir::Right => {
            // Right of the deletion must be a match state
            let Some(StateRange::M(right1)) = right1 else {
                return IdxAdjustment::next();
            };

            // The match state must have sufficient length to cover the new
            // `to` range
            let distance = rule.distance();
            // We require a strictly greater match length to prevent
            // pathological cases such as shifting to edge of alignment,
            // shifting adjacent to another deletion, etc. TODO: Convert to
            // strict inequality.
            if right1.len() <= distance {
                return IdxAdjustment::next();
            }

            *any_rewritten = true;

            // Shrink the match state and shift the deletion
            right1.cut_start(distance);
            del.shift_right(distance);

            if let Some(StateRange::M(left1)) = left1 {
                // Extend the left match state
                left1.extend_end(distance);

                IdxAdjustment::next()
            } else {
                // Validity: this will be in bounds for the query and
                // reference based on the cut_start and shift_right calls
                // above (this is effectively undoing that on one side)
                let query_range = right1.query_range.start - distance..right1.query_range.start;
                let ref_range = del.ref_range.start - distance..del.ref_range.start;

                let new_match_state = StateRange::M(MatchRange { query_range, ref_range });

                // Insert before the deletion
                IdxAdjustment {
                    advance: true,
                    edits:   Some(StateVecEdits {
                        insert_left: Some(new_match_state),
                        ..Default::default()
                    }),
                }
            }
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
}
