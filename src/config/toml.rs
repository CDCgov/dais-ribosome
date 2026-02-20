use serde::Deserialize;
use serde_derive::Deserialize;
use std::{
    collections::HashMap,
    path::{Path, PathBuf},
};
use zoe::data::err::ResultWithErrorContext;

/// Root configuration structure parsed from `modules.toml`.
#[derive(Debug, Clone, Deserialize)]
pub struct TomlConfig {
    #[serde(rename = "module")]
    pub modules: Vec<ConfiguredModule>,
}

impl TomlConfig {
    /// Load configuration from a TOML file.
    pub fn from_file(path: &Path) -> std::io::Result<Self> {
        let raw_toml = std::fs::read_to_string(path)?;
        Ok(toml::from_str::<TomlConfig>(&raw_toml).with_type_context::<TomlConfig>()?)
    }

    /// Find a module by name, returning it along with the other modules'
    /// metadata.
    ///
    /// This is used when forming [`ModuleData`]. The other modules' metadata
    /// are used to populate `other_modules`, which is used in warning messages
    /// in [`print_unimplemented_ctypes`].
    pub fn find_module(self, name: &str, module_path: &Path) -> Option<(ConfiguredModule, Vec<(String, PathBuf)>)> {
        let mut selected = None;
        let mut others = Vec::new();

        for m in self.modules {
            if m.name == name {
                selected = Some(m);
            } else {
                let ref_path = module_path.join(&m.name).join(&m.references);
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
    pub alignment:  HashMap<String, AlignmentParams>,
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

/// Alignment scoring parameters.
#[derive(Debug, Clone, Deserialize)]
pub struct AlignmentParams {
    #[serde(rename = "match")]
    pub match_score: i8,
    pub mismatch:    i8,
    #[serde(deserialize_with = "deserialize_gap_penalty")]
    pub gap_open:    i8,
    #[serde(deserialize_with = "deserialize_gap_penalty")]
    pub gap_extend:  i8,
}

// TODO: Currently it is not validating the range...

/// Deserializes a gap penalty, validating the range and normalizing to
/// negative.
fn deserialize_gap_penalty<'de, D>(deserializer: D) -> Result<i8, D::Error>
where
    D: serde::Deserializer<'de>, {
    let value: i8 = Deserialize::deserialize(deserializer)?;
    Ok(-value.abs())
}
