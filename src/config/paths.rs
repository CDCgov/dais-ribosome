use std::{io::Error, path::PathBuf};

const MODULE_TOML_SUFFIX: &str = "ribosome_res/modules.toml";

/// Locate the `modules.toml` configuration file.
pub fn find_modules_toml() -> Result<PathBuf, Error> {
    let exe_path = std::env::current_exe()?;
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
