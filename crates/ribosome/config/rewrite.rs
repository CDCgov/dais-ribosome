//! The data structures and parsing for the deletion repositioning TOML file.

use crate::{
    data::keys::RefKey,
    ranges::{MatchRange, ParseOneBasedInclusive, StateRange},
};
use serde::{Deserialize, Deserializer, de::Error};
use std::{collections::HashMap, ops::Range, path::Path};
use zoe::data::err::ResultWithErrorContext;

/// The rewrite rules for a module, indexed by reference ID and compound type.
#[derive(Debug, Default)]
pub struct RewriteRules {
    /// The rules for rewriting a deletion.
    pub(crate) deletions: HashMap<RefKey, Vec<RewriteRanges>>,
}

/// The ranges for a rewrite rule.
///
/// The `from` and `to` range will be the same length and will be not equal.
#[derive(Debug)]
pub(crate) struct RewriteRanges {
    /// The original range where the indel must be present to apply the rule, in
    /// 0-based reference nucleotide coordinates.
    from: Range<usize>,
    /// The destination range where the deletion will be positioned.
    to:   Range<usize>,
}

// TODO: support a `[cds]` section alongside `[genome]`.
/// A helper type for deserializing [`RewriteRanges`].
#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct RewriteRulesRaw {
    #[serde(default)]
    genome: GenomeRulesRaw,
}

// TODO: support `[[genome.insertions]]`.
#[derive(Debug, Default, Deserialize)]
#[serde(deny_unknown_fields)]
struct GenomeRulesRaw {
    #[serde(default)]
    deletions: Vec<DeletionEntry>,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct DeletionEntry {
    rule: DeletionRule,
}

#[derive(Debug)]
struct DeletionRule {
    ref_key: RefKey,
    ranges:  RewriteRanges,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct DeletionRuleRaw {
    ctype:        String,
    reference_id: String,
    from:         String,
    to:           String,
}

impl<'de> Deserialize<'de> for DeletionRule {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>, {
        let DeletionRuleRaw {
            ctype,
            reference_id,
            from: from_str,
            to: to_str,
        } = DeletionRuleRaw::deserialize(deserializer)?;

        let Ok(from) = Range::<usize>::parse_inclusive(&from_str) else {
            return Err(D::Error::custom(format!(
                "failed to parse range {from_str}. Expected <start>..<end>"
            )));
        };

        let Ok(to) = Range::<usize>::parse_inclusive(&to_str) else {
            return Err(D::Error::custom(format!(
                "failed to parse range {to_str}. Expected <start>..<end>"
            )));
        };

        if from.len() != to.len() {
            return Err(D::Error::custom(format!(
                "the from ({from_str}) and to ({to_str}) ranges are not the same length"
            )));
        }

        if from == to {
            return Err(D::Error::custom(format!(
                "the same range was found for both from and to: {from_str}"
            )));
        }

        // TODO: Validate no tabs
        // Validity: we confirmed the ranges are the same length and not equal
        Ok(DeletionRule {
            ref_key: RefKey::new(reference_id, ctype),
            ranges:  RewriteRanges { from, to },
        })
    }
}

impl RewriteRules {
    /// Reads the rewrite rules TOML file from the specified `path`.
    ///
    /// ## Errors
    ///
    /// IO and parsing errors are propagated with path context.
    pub fn from_file(path: impl AsRef<Path>) -> std::io::Result<Self> {
        let path = path.as_ref();

        let raw = std::fs::read_to_string(path)
            .and_then(|raw_toml| Ok(toml::from_str::<RewriteRulesRaw>(&raw_toml).with_type_context::<RewriteRulesRaw>()?))
            .with_path_context("Failed to parse rewrite rules TOML file", path)?;

        let mut deletions: HashMap<RefKey, Vec<RewriteRanges>> = HashMap::new();

        for DeletionEntry { rule } in raw.genome.deletions {
            deletions.entry(rule.ref_key).or_default().push(rule.ranges);
        }

        Ok(RewriteRules { deletions })
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
    pub(crate) fn rewrite_deletion(&self, genome_aln_states: &mut Vec<StateRange>) -> bool {
        // Find a deletion whose range equals `from`
        let Some(del_idx) = genome_aln_states
            .iter()
            .position(|state| matches!(state, StateRange::D(del) if del.ref_range == self.from))
        else {
            return false;
        };

        // Get the flanking states
        let Some(SurroundingStates {
            left,
            current: StateRange::D(del),
            right,
        }) = get_surrounding_states(del_idx, genome_aln_states)
        else {
            // This should be unreachable, since del_idx is known to be in
            // bounds and correspond to a deletion
            return false;
        };

        match self.direction() {
            RepositionDir::Left => {
                // Left of the deletion must be a match state
                let Some(StateRange::M(left)) = left else { return false };

                // The match state must have sufficient length to cover the new
                // `to` range
                let distance = self.distance();
                // We require a strictly greater match length to prevent
                // pathological cases such as shifting to edge of alignment,
                // shifting adjacent to another deletion, etc. TODO: Convert to
                // strict inequality.
                if left.len() <= distance {
                    return false;
                }

                // Shrink the match state and shift the deletion
                left.cut_end(distance);
                del.shift_left(distance);

                if let Some(StateRange::M(right)) = right {
                    // Extend the right match state
                    right.extend_start(distance);
                } else {
                    // Validity: this will be in bounds for the query and
                    // reference based on the cut_end and shift_left calls above
                    // (this is effectively undoing that on one side)
                    let query_range = left.query_range.end..left.query_range.end + distance;
                    let ref_range = del.ref_range.end..del.ref_range.end + distance;

                    // Insert after the deletion
                    genome_aln_states.insert(del_idx + 1, StateRange::M(MatchRange { query_range, ref_range }));
                }
            }
            RepositionDir::Right => {
                // Right of the deletion must be a match state
                let Some(StateRange::M(right)) = right else { return false };

                // The match state must have sufficient length to cover the new
                // `to` range
                let distance = self.distance();
                // We require a strictly greater match length to prevent
                // pathological cases such as shifting to edge of alignment,
                // shifting adjacent to another deletion, etc. TODO: Convert to
                // strict inequality.
                if right.len() <= distance {
                    return false;
                }

                // Shrink the match state and shift the deletion
                right.cut_start(distance);
                del.shift_right(distance);

                if let Some(StateRange::M(left)) = left {
                    // Extend the left match state
                    left.extend_end(distance);
                } else {
                    // Validity: this will be in bounds for the query and
                    // reference based on the cut_start and shift_right calls
                    // above (this is effectively undoing that on one side)
                    let query_range = right.query_range.start - distance..right.query_range.start;
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

/// A helper struct used in [`rewrite_deletion`], holding simultaneous mutable
/// references to a state and its flanking states if present.
///
/// [`rewrite_deletion`]: RewriteRanges::rewrite_deletion
struct SurroundingStates<'a> {
    /// The state left of the specified index, if one exists.
    left:    Option<&'a mut StateRange>,
    /// The state at the specified index.
    current: &'a mut StateRange,
    /// The state right of the specified index, if one exists.
    right:   Option<&'a mut StateRange>,
}

/// A helper function to extract simultaneous mutable references to a state at
/// `idx` as well as its flanking states.
///
/// `None` is returned if `idx` is out of bounds.
fn get_surrounding_states(idx: usize, ranges: &mut [StateRange]) -> Option<SurroundingStates<'_>> {
    let (left, current_and_right) = ranges.split_at_mut_checked(idx)?;
    let (current, right) = current_and_right.split_first_mut()?;

    let left = left.last_mut();
    let right = right.first_mut();

    Some(SurroundingStates { left, current, right })
}
