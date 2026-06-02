//! The data structures and parsing for the TOML configuration file.

use serde::{
    Deserialize, Deserializer,
    de::{Error, Unexpected},
};
use std::{
    collections::{HashMap, hash_map::Entry},
    path::{Path, PathBuf},
};
use zoe::data::{WeightMatrix, err::ResultWithErrorContext};

/// Root configuration structure parsed from `modules.toml`.
///
/// ## Notes
///
/// The alignment parameters specified for mismatch, gap open, and gap extend
/// are automatically converted to non-positive.
#[derive(Clone, Debug)]
pub struct TomlConfig {
    /// The available modules.
    ///
    /// All modules should have distinct names, ignoring case. This invariant is
    /// upheld by the [`Deserialize`] implementation.
    pub modules: Vec<ConfiguredModule>,
}

impl TomlConfig {
    /// Load configuration from a TOML file.
    ///
    /// ## Errors
    ///
    /// If any IO errors or parsing errors occur, an error with path context
    /// (and type context for parsing errors) is returned.
    pub fn from_file(path: impl AsRef<Path>) -> std::io::Result<Self> {
        let path = path.as_ref();

        Ok(std::fs::read_to_string(path)
            .and_then(|raw_toml| Ok(toml::from_str::<TomlConfig>(&raw_toml).with_type_context::<TomlConfig>()?))
            .with_path_context("Failed to parse TOML file", path)?)
    }

    /// Find a module by name, returning it along with the paths to other
    /// modules' references.
    ///
    /// This is used when forming [`AnnotationModule`]. The second return value
    /// is a vector containing tuples with the module names and the paths to the
    /// reference sequences. These are used to populate `other_modules`.
    ///
    /// [`AnnotationModule`]: crate::config::annotation_module::AnnotationModule
    pub(crate) fn find_module(
        &self, name: &str, modules_dir: &Path,
    ) -> Option<(&ConfiguredModule, Vec<(&String, PathBuf)>)> {
        let mut selected = None;
        let mut others = Vec::new();

        for m in &self.modules {
            if m.name.eq_ignore_ascii_case(name)
                || m.alternative_names.iter().any(|alt_name| name.eq_ignore_ascii_case(alt_name))
            {
                selected = Some(m);
            } else {
                let ref_path = modules_dir.join(&m.name).join(&m.references);
                others.push((&m.name, ref_path));
            }
        }

        selected.map(|m| (m, others))
    }
}

/// A helper type for parsing [`TomlConfig`], allowing us to validate that the
/// names of the modules are unique.
#[derive(Deserialize)]
struct TomlConfigRaw {
    #[serde(rename = "module")]
    modules: Vec<ConfiguredModule>,
}

impl<'de> Deserialize<'de> for TomlConfig {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>, {
        let raw = TomlConfigRaw::deserialize(deserializer)?;

        // The currently-encountered names as keys. If None, it is a module
        // name. If Some, it is an alternative name for the specified module
        let mut names: HashMap<String, Option<String>> = HashMap::new();

        for module in &raw.modules {
            let name = module.name.to_ascii_lowercase();

            match names.entry(name.clone()) {
                Entry::Occupied(entry) => {
                    return match entry.get() {
                        Some(other_module_name) => Err(D::Error::custom(format!(
                            "{name} appeared as both a module and as an alternative name for {other_module_name}"
                        ))),
                        None => Err(D::Error::custom(format!("two modules with name {name} were found"))),
                    };
                }
                Entry::Vacant(entry) => {
                    entry.insert(None);
                }
            }

            for alt_name in &module.alternative_names {
                let alt_name = alt_name.to_ascii_lowercase();

                match names.entry(alt_name.clone()) {
                    Entry::Occupied(entry) => {
                        return match entry.get() {
                            Some(other_module_name) => Err(D::Error::custom(format!(
                                "{alt_name} appeared as an alternative name for both {name} and {other_module_name}"
                            ))),
                            None => Err(D::Error::custom(format!(
                                "{alt_name} exists as a module and as an alternative name for {name}"
                            ))),
                        };
                    }
                    Entry::Vacant(entry) => {
                        entry.insert(Some(module.name.clone()));
                    }
                }
            }
        }

        Ok(TomlConfig { modules: raw.modules })
    }
}

/// Configuration for a single annotation module (e.g., `flu`, `cov`, or `rsv`).
#[derive(Clone, Debug, Deserialize)]
pub struct ConfiguredModule {
    /// The name of the module (e.g., `flu`, `cov`, or `rsv`). This must
    /// correspond to a folder in `ribosome_res`.
    pub name:              String,
    /// An optional version for the module (e.g., `2.0-alpha`).
    pub version:           Option<String>,
    /// Any alternative names that can be used to refer to the module.
    #[serde(default)]
    pub alternative_names: Vec<String>,
    /// The file name for the FASTA file containing the references. This should
    /// be a relative path within the module folder.
    pub references:        PathBuf,
    // The path containing the codon position weights.
    pub weights:           PathBuf,
    /// The file name for the TSV file containing the coding sequence (CDS)
    /// specifications. This should be a relative path within the module folder.
    pub cds_spec:          PathBuf,
    /// The alignment method to use.
    #[serde(default)]
    pub alignment_method:  AlignmentMethod,
    /// The output formatting options.
    pub formatting:        Formatting,
    /// Rules allowing customization of the annotation process.
    pub rules:             Rules,
    /// Collection of alignment weights for the module (and specific compound
    /// types within it).
    pub alignment:         AlignmentWeights,
}

/// The supported alignment methods in DAIS-ribosome.
#[derive(Copy, Clone, Eq, PartialEq, Debug, Default, Deserialize)]
pub enum AlignmentMethod {
    /// The 1-pass algorithm, which is faster for smaller sequences (e.g., flu).
    #[serde(rename = "one-pass")]
    OnePass,
    /// The 3-pass algorithm, which is more memory efficient for larger
    /// sequences (e.g., covid).
    #[default]
    #[serde(rename = "three-pass")]
    ThreePass,
}

/// Collection of alignment weights for a module.
///
/// Provides default scoring parameters plus optional overrides for specific
/// compound types.
#[derive(Clone, Debug, Deserialize)]
pub struct AlignmentWeights {
    /// Default scoring matrix and gap parameters.
    pub default: AlignmentParams,

    /// Per-compound-type overrides.
    #[serde(flatten)]
    pub overrides: HashMap<String, AlignmentParams>,
}

impl AlignmentWeights {
    /// Gets the scoring parameters for a compound type.
    ///
    /// Returns the override for the compound type if it exists, otherwise
    /// returns defaults.
    pub fn get(&self, compound_type: &str) -> &AlignmentParams {
        self.overrides.get(compound_type).unwrap_or(&self.default)
    }
}

/// Output formatting options for a module.
///
/// Currently, all of these are related to padding on the right, which is not
/// required for ensuring the proper reading frame, but may be useful to
/// downstream applications which expect sequences of a given length (or the
/// same length, such as multiple sequence alignment programs).
#[derive(Debug, Clone, Deserialize)]
pub struct Formatting {
    /// Whether to add right padding to the amino acid sequence.
    #[serde(default = "pad_default")]
    pub right_pad_aa:  bool,
    /// Whether to add right padding to the coding sequence.
    #[serde(default = "pad_default")]
    pub right_pad_cds: bool,
    /// Whether to add right padding to the genome alignment.
    #[serde(default = "pad_default")]
    pub right_pad_gen: bool,
}

/// Returns the default value for whether to perform padding. This is used by
/// `serde` in [`Formatting`].
const fn pad_default() -> bool {
    true
}

/// Rules allowing customization of the annotation process.
#[derive(Debug, Clone, Deserialize)]
pub struct Rules {
    /// If the genome alignment reaches the end of the reference but does not
    /// end in a stop codon as expected, then this rule causes the alignment to
    /// be extended until the next in-frame stop codon. The extra bases are
    /// represented as an insertion at the end of the genome and any products
    /// extending to the end of the genome.
    #[serde(default)]
    pub list_contig_stop_extension: bool,

    /// If the query does not start with a start codon, this rule will trim the
    /// start of the query up to this start codon before carrying out alignment.
    /// However, this trimming is only applied if the length of the query is at
    /// least the length of the reference (i.e., it is feasible that the query
    /// just has extra non-coding state on the left compared to the reference).
    ///
    /// This rule is designed to help reduce the complexity of the alignment
    /// that is formed. The trimmed bases may be re-added if the
    /// `repairable_end_limit` rule is enabled.
    #[serde(default)]
    pub chew_to_start: bool,

    /// If there are soft clipped nucleotides within the specified limit on
    /// either side of the alignment, then extend the alignment with match
    /// states to add them. Setting this to 0 disables the rule.
    ///
    /// When mismatches are present near either end of the query/reference,
    /// local alignment can cause the ends of the query/reference to not be
    /// included. This method _decreases_ the optimality of the alignment from a
    /// Smith-Waterman standpoint, but may produce better products or genome
    /// alignments.
    ///
    /// Specifically, to extend the alignment on a side, either the number of
    /// clipped bases in the query _or_ the number of clipped bases in the
    /// reference must be at most the specified limit.
    ///
    /// This rule may add back the bases removed by `chew_to_start`.
    #[serde(default)]
    pub repairable_end_limit: usize,
}

/// The alignment scoring parameters.
///
/// ## Validity
///
/// At least one of `gap_open` or `gap_extend` must be nonzero, to ensure that
/// optimal local alignments do not begin or end with indels.
#[derive(Debug, Clone)]
pub struct AlignmentParams {
    /// The weight matrix for alignment.
    ///
    /// ## Validity
    ///
    /// The mismatch score must be non-positive, and cannot be -128.
    pub matrix:     WeightMatrix<'static, i8, 5>,
    /// The score for opening a gap.
    ///
    /// ## Validity
    ///
    /// This must be non-positive, and cannot be -128.
    pub gap_open:   i8,
    /// The score for extending a gap.
    ///
    /// ## Validity
    ///
    /// This must be non-positive, and cannot be -128.
    pub gap_extend: i8,
}

/// A helper type for parsing [`AlignmentParams`] that does not impose any
/// conditions or checking on the integers.
#[derive(Deserialize)]
struct AlignmentParamsRaw {
    #[serde(rename = "match")]
    match_score: i8,

    #[serde(deserialize_with = "deserialize_penalty")]
    mismatch: i8,

    #[serde(deserialize_with = "deserialize_penalty")]
    gap_open: i8,

    #[serde(deserialize_with = "deserialize_penalty")]
    gap_extend: i8,
}

impl<'de> Deserialize<'de> for AlignmentParams {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>, {
        // Validity: AlignmentParamsRaw handles converting penalties to
        // non-positive and ensuring they are not -128.
        let raw = AlignmentParamsRaw::deserialize(deserializer)?;

        // Validity: Verify validity requirement on AlignmentParams
        if raw.gap_open == 0 && raw.gap_extend == 0 {
            return Err(D::Error::custom("gap_open and gap_extend cannot both be 0"));
        }

        let matrix = WeightMatrix::new_dna_matrix(raw.match_score, raw.mismatch, Some(b'N'));

        Ok(AlignmentParams {
            matrix,
            gap_open: raw.gap_open,
            gap_extend: raw.gap_extend,
        })
    }
}

/// Deserializes an alignment penalty, validating the range and normalizing to
/// be non-positive.
///
/// ## Errors
///
/// If the penalty is -128, then an error is thrown since this is out of the
/// range *Zoe* can handle.
fn deserialize_penalty<'de, D>(deserializer: D) -> Result<i8, D::Error>
where
    D: serde::Deserializer<'de>, {
    let value: i8 = Deserialize::deserialize(deserializer)?;
    match value {
        -128 => Err(D::Error::invalid_value(
            Unexpected::Signed(value as i64),
            &"an integer of absolute value at most 127",
        )),
        -127..0 => Ok(value),
        0.. => Ok(-value),
    }
}
