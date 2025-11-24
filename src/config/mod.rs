//! Configuration loading and path resolution.
//!
//! This module handles:
//! - Locating the `modules.toml` configuration file
//! - Parsing module configurations from TOML
//! - Resolving paths to module resources

mod paths;
mod toml;

pub use paths::{find_modules_toml, module_resource_dir};
pub use toml::{AlignmentParams, ConfiguredModule, Formatting, Rules, TomlConfig, suggest_module_for_compound_type};
