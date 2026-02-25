use serde::{
    Deserialize, Deserializer,
    de::{Error, Unexpected},
};
use serde_derive::Deserialize;
use std::{
    collections::HashMap,
    path::{Path, PathBuf},
};
use zoe::data::{WeightMatrix, err::ResultWithErrorContext};

/// Root configuration structure parsed from `modules.toml`.
#[derive(Debug, Clone, Deserialize)]
pub struct TomlConfig {
    #[serde(rename = "module")]
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

    /// Find a module by name, returning it along with the other modules'
    /// metadata.
    ///
    /// This is used when forming [`ModuleData`]. The second return value is a
    /// vector containing tuples with the module names and the paths to the
    /// reference sequences. These are used to populate `other_modules`, which
    /// is used in warning messages in [`print_unimplemented_ctypes`].
    pub fn find_module(self, name: &str, modules_dir: &Path) -> Option<(ConfiguredModule, Vec<(String, PathBuf)>)> {
        let mut selected = None;
        let mut others = Vec::new();

        for m in self.modules {
            if m.name == name {
                selected = Some(m);
            } else {
                let ref_path = modules_dir.join(&m.name).join(&m.references);
                others.push((m.name, ref_path));
            }
        }

        selected.map(|m| (m, others))
    }
}

/// Configuration for a single annotation module (e.g., `flu`, `cov`, or `rsv`).
#[derive(Clone, Debug, Deserialize)]
pub struct ConfiguredModule {
    /// The name of the module (e.g., `flu`, `cov`, or `rsv`). This must
    /// correspond to a folder in `ribosome_res`.
    pub name:       String,
    /// An optional version for the module (e.g., `2.0-alpha`).
    pub version:    Option<String>,
    /// The file name for the FASTA file containing the references. This should
    /// be a relative path within the module folder.
    pub references: PathBuf,
    // TODO: What is this?
    pub weights:    PathBuf,
    /// The file name for the TSV file containing the coding sequence (CDS)
    /// specifications. This should be a relative path within the module folder.
    pub cds_spec:   PathBuf,
    // TODO: What is this?
    pub formatting: Formatting,
    // TODO: What is this?
    pub rules:      Rules,
    pub alignment:  AlignmentWeights,
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

// TODO: What are these fields?

/// Output formatting options for a module.
#[derive(Debug, Clone, Deserialize)]
pub struct Formatting {
    #[serde(default = "pad_default")]
    pub right_pad_aa:  bool,
    #[serde(default = "pad_default")]
    pub right_pad_cds: bool,
    #[serde(default = "pad_default")]
    pub right_pad_gen: bool,
}

/// Returns the default value for whether to perform padding. This is used by
/// `serde` in [`Formatting`].
const fn pad_default() -> bool {
    true
}

#[derive(Debug, Clone, Deserialize)]
pub struct Rules {
    #[serde(default)]
    pub list_contig_stop_extension: bool,
    #[serde(default)]
    pub chew_to_start:              bool,
    pub repairable_end_limit:       Option<usize>,
}

// TODO: Should we be negating the mismatch penalty too when parsing?

/// A helper type for parsing [`AlignmentParams`] that does not impose any
/// conditions or checking on the integers.
#[derive(Deserialize)]
struct AlignmentParamsRaw {
    #[serde(rename = "match")]
    match_score: i8,
    mismatch:    i8,
    #[serde(deserialize_with = "deserialize_gap_penalty")]
    gap_open:    i8,
    #[serde(deserialize_with = "deserialize_gap_penalty")]
    gap_extend:  i8,
}

/// Alignment scoring parameters.
#[derive(Debug, Clone)]
pub struct AlignmentParams {
    pub matrix:     WeightMatrix<'static, i8, 5>,
    pub gap_open:   i8,
    pub gap_extend: i8,
}

impl<'de> Deserialize<'de> for AlignmentParams {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>, {
        let raw = AlignmentParamsRaw::deserialize(deserializer)?;

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

/// Deserializes a gap penalty, validating the range and normalizing to
/// negative.
///
/// ## Errors
///
/// If the gap penalty is -128, then an error is thrown since this is out of the
/// range *Zoe* can handle.
fn deserialize_gap_penalty<'de, D>(deserializer: D) -> Result<i8, D::Error>
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
