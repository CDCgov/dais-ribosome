//! Configuration loading and path resolution.
//!
//! This module handles:
//! - Locating the `modules.toml` configuration file
//! - Parsing module configurations from TOML
//! - Resolving paths to module resources

mod toml;
pub use toml::*;

use std::{io::Error, path::PathBuf};
use zoe::data::err::ResultWithErrorContext;

const MODULE_TOML_SUFFIX: &str = "ribosome_res/modules.toml";

/// Similar to [`std::env::current_exe`], but includes context for the error.
pub fn current_exe() -> std::io::Result<PathBuf> {
    Ok(std::env::current_exe().with_context("Failed to find the path to the executable")?)
}

/// Locates the `modules.toml` configuration file, returning a full filesystem
/// path.
///
/// This must be located directly within the `ribosome_res` folder. The
/// `ribosome_res` folder is searched for in the same folder as the current
/// executable, as well as the grandparent folder (to support `target/debug` or
/// `target/release` builds).
pub fn find_modules_toml() -> std::io::Result<PathBuf> {
    let exe_path = current_exe()?;
    let Some(exe_dir) = exe_path.parent() else {
        return Err(Error::other(
            "The executable should have a parent folder, but somehow did not.",
        ));
    };

    // Check same directory as executable
    let candidate = exe_dir.join(MODULE_TOML_SUFFIX);
    if candidate.exists() {
        return Ok(candidate);
    }

    // Check grandparent (for target/debug or target/release builds)
    if let Some(grandparent) = exe_dir.parent().and_then(|p| p.parent()) {
        let candidate = grandparent.join(MODULE_TOML_SUFFIX);
        if candidate.exists() {
            return Ok(candidate);
        }
    }

    Err(Error::other("Could not find 'modules.toml' file"))
}
