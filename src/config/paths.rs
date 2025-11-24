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

/// Returns the directory containing a module's resource files.
pub fn module_resource_dir(modules_toml_path: &std::path::Path, module_name: &str) -> PathBuf {
    modules_toml_path
        .parent()
        .map(std::path::Path::to_path_buf)
        .unwrap_or_else(|| PathBuf::from("."))
        .join(module_name)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_module_resource_dir() {
        let toml_path = PathBuf::from("/app/ribosome_res/modules.toml");
        let dir = module_resource_dir(&toml_path, "flu");
        assert_eq!(dir, PathBuf::from("/app/ribosome_res/flu"));
    }
}
