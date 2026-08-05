//! The data structures and parsing for the deletion repositioning TOML file.

use crate::{data::keys::RefKey, ranges::ParseOneBasedInclusive};
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
#[derive(Eq, PartialEq, Debug)]
pub(crate) struct RewriteRanges {
    /// The original range where the indel must be present to apply the rule, in
    /// 0-based reference nucleotide coordinates.
    pub from: Range<usize>,
    /// The destination range where the deletion will be positioned.
    pub to:   Range<usize>,
}

impl Ord for RewriteRanges {
    /// This method returns an [`Ordering`] between `self` and `other`.
    ///
    /// The sort order is based first on the starting index of `from`, then the
    /// ending index of `from`, the starting index of `to`, and the ending index
    /// of `to`.
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.from
            .start
            .cmp(&other.from.start)
            .then_with(|| self.from.end.cmp(&other.from.end))
            .then_with(|| self.to.start.cmp(&other.to.start))
            .then_with(|| self.to.end.cmp(&other.to.end))
    }
}

impl PartialOrd for RewriteRanges {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
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
