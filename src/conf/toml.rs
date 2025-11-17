use serde::Deserialize;
use serde_derive::Deserialize;
use std::{
    collections::HashMap,
    path::{Path, PathBuf},
};

pub fn process_toml(path: &Path) -> TomlConfig {
    let raw_toml = std::fs::read_to_string(path).unwrap();
    let config = toml::from_str::<TomlConfig>(&raw_toml).unwrap();
    eprintln!("{modules:#?}", modules = config.modules);
    config
}

#[derive(Debug, Deserialize)]
pub struct TomlConfig {
    #[serde(rename = "module")]
    pub modules: Vec<Module>,
}

#[derive(Debug, Deserialize)]
pub struct Module {
    pub name:       String,
    pub version:    Option<String>,
    pub references: PathBuf,
    pub weights:    PathBuf,
    pub cds_spec:   PathBuf,
    pub formatting: Formatting,
    pub alignment:  HashMap<String, AlignmentParams>,
}

#[derive(Debug, Deserialize)]
pub struct Formatting {
    right_pad_aa:  bool,
    right_pad_cds: bool,
}

#[derive(Debug, Deserialize)]
pub struct AlignmentParams {
    #[serde(rename = "match")]
    pub match_score: i8,
    pub mismatch:    i8,
    pub gap_open:    i8,
    pub gap_extend:  i8,
}

#[derive(Debug, Deserialize)]
pub struct Products {
    pub protein: String,
    #[serde(deserialize_with = "deserialize_ranges")]
    pub ranges:  Vec<std::ops::Range<usize>>,
}

fn deserialize_ranges<'de, D>(deserializer: D) -> Result<Vec<std::ops::Range<usize>>, D::Error>
where
    D: serde::Deserializer<'de>, {
    let raw_ranges: Vec<(usize, usize)> = Deserialize::deserialize(deserializer)?;
    Ok(raw_ranges.into_iter().map(|(start, end)| (start - 1)..end).collect())
}
