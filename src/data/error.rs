//! Module loading errors.

use std::{fmt, path::PathBuf};

/// Error type for profile building failures.
pub type ProfileBuildError = zoe::alignment::ProfileError;

/// Errors that can occur when loading a module.
#[derive(Debug)]
pub enum ModuleLoadError {
    /// Failed to read a file.
    Io { path: PathBuf, source: std::io::Error },
    /// Requested module not found in configuration.
    ModuleNotFound { name: String },
    /// Module configuration is invalid.
    InvalidConfig { module: String, reason: &'static str },
    /// Failed to build alignment profiles.
    ProfileBuild { source: ProfileBuildError },
    /// Could not locate modules.toml configuration file.
    ConfigNotFound,
    /// Data validation error.
    Validation { message: String },
}

impl ModuleLoadError {
    /// Create an IO error with path context.
    pub fn io(path: &std::path::Path, source: std::io::Error) -> Self {
        ModuleLoadError::Io {
            path: path.to_path_buf(),
            source,
        }
    }

    /// Create a module not found error.
    pub fn module_not_found(name: impl Into<String>) -> Self {
        ModuleLoadError::ModuleNotFound { name: name.into() }
    }

    /// Create an invalid configuration error.
    pub fn invalid_config(module: impl Into<String>, reason: &'static str) -> Self {
        ModuleLoadError::InvalidConfig {
            module: module.into(),
            reason,
        }
    }

    /// Create a validation error.
    pub fn validation(message: impl Into<String>) -> Self {
        ModuleLoadError::Validation { message: message.into() }
    }
}

impl fmt::Display for ModuleLoadError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ModuleLoadError::Io { path, source: _ } => {
                write!(f, "IO failures in file: {}", path.display())
            }
            ModuleLoadError::ModuleNotFound { name } => {
                write!(f, "Module '{name}' not found in configuration")
            }
            ModuleLoadError::InvalidConfig { module, reason } => {
                write!(f, "Module '{module}' has invalid configuration: {reason}")
            }
            ModuleLoadError::ProfileBuild { source: _ } => {
                write!(f, "Failed to build alignment profiles")
            }
            ModuleLoadError::ConfigNotFound => {
                write!(f, "Could not locate modules.toml configuration file")
            }
            ModuleLoadError::Validation { message } => {
                write!(f, "Validation error: {message}")
            }
        }
    }
}

impl std::error::Error for ModuleLoadError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            ModuleLoadError::Io { source, .. } => Some(source),
            ModuleLoadError::ProfileBuild { source } => Some(source),
            _ => None,
        }
    }
}

impl From<ProfileBuildError> for ModuleLoadError {
    fn from(source: ProfileBuildError) -> Self {
        ModuleLoadError::ProfileBuild { source }
    }
}

impl zoe::data::err::GetCode for ModuleLoadError {
    fn get_code(&self) -> i32 {
        if let ModuleLoadError::Io { source, .. } = self {
            source.get_code()
        } else {
            1
        }
    }
}
