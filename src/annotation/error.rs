//! Error construction helpers for data loading.

use zoe::data::err::{ErrorWithContext, GetCode};

/// A minimal enum containing error variants which must be matched on in
/// DAIS-ribosome. All other errors should use [`RibosomeError::Io`].
#[derive(Debug)]
pub enum RibosomeError {
    UnimplementedCtype(String),
    EmptyFile(std::path::PathBuf),
    Io(std::io::Error),
}

impl std::error::Error for RibosomeError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            RibosomeError::Io(e) => e.source(),
            _ => None,
        }
    }
}

impl GetCode for RibosomeError {
    fn get_code(&self) -> i32 {
        match self {
            RibosomeError::Io(e) => e.get_code(),
            RibosomeError::EmptyFile(_) => 66, // EX_NOINPUT
            _ => 1,
        }
    }
}

impl std::fmt::Display for RibosomeError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            RibosomeError::UnimplementedCtype(c) => write!(f, "Unimpelmented ctype '{c}'"),
            RibosomeError::EmptyFile(p) => write!(f, "Empty file: {}", p.display()),
            RibosomeError::Io(e) => write!(f, "{e}"),
        }
    }
}

impl From<&str> for RibosomeError {
    fn from(value: &str) -> Self {
        RibosomeError::Io(std::io::Error::other(value))
    }
}

impl From<String> for RibosomeError {
    fn from(value: String) -> Self {
        RibosomeError::Io(std::io::Error::other(value))
    }
}

impl From<std::io::Error> for RibosomeError {
    fn from(e: std::io::Error) -> Self {
        RibosomeError::Io(e)
    }
}

impl From<ErrorWithContext> for RibosomeError {
    fn from(e: ErrorWithContext) -> Self {
        RibosomeError::Io(std::io::Error::other(e))
    }
}
