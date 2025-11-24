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

    /// Find a module by name, returning it along with the other modules' metadata.
    pub fn find_module(self, name: &str) -> Option<(ConfiguredModule, Vec<(String, PathBuf)>)> {
        let mut selected = None;
        let mut others = Vec::new();

        for m in self.modules {
            if m.name == name {
                selected = Some(m);
            } else {
                others.push((m.name, m.references));
            }
        }

        selected.map(|m| (m, others))
    }
}

/// Configuration for a single annotation module (e.g., flu, cov, rsv).
#[derive(Debug, Clone, Deserialize)]
pub struct ConfiguredModule {
    pub name:       String,
    pub version:    Option<String>,
    pub references: PathBuf,
    pub weights:    PathBuf,
    pub cds_spec:   PathBuf,
    pub formatting: Formatting,
    pub rules:      Rules,
    pub alignment:  HashMap<String, AlignmentParams>,
}

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

/// Deserializes a gap penalty, validating the range and normalizing to negative.
fn deserialize_gap_penalty<'de, D>(deserializer: D) -> Result<i8, D::Error>
where
    D: serde::Deserializer<'de>, {
    let value: i8 = Deserialize::deserialize(deserializer)?;
    match value {
        0 => Ok(0),
        1.. => Ok(-value),
        ..=-1 => Ok(value),
    }
}

/// Suggest which module might contain a given compound type.
pub fn suggest_module_for_compound_type<'a>(
    modules: &'a [ConfiguredModule], compound_type: &str, exclude_module: &str,
) -> Option<&'a str> {
    modules
        .iter()
        .filter(|module| module.name != exclude_module)
        .find(|module| module.alignment.contains_key(compound_type))
        .map(|module| module.name.as_str())
}
