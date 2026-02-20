//! Configuration loading and path resolution.
//!
//! This module handles:
//! - Locating the `modules.toml` configuration file
//! - Parsing module configurations from TOML
//! - Resolving paths to module resources

mod paths;
mod toml;

pub use paths::{current_exe, find_modules_toml};
pub use toml::{AlignmentParams, ConfiguredModule, Formatting, Rules, TomlConfig};
