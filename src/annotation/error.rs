//! Error construction helpers for data loading.

use zoe::data::err::{ErrorWithContext, GetCode};

#[derive(Debug)]
pub enum RibosomeError {
    InvalidFastaFormat,
    InvalidTsvFormat,
    BlankFirstLine(std::path::PathBuf),
    InvalidSequence(String),
    NoCtype(String),
    UnimplementedCtype(String),
    /// Query could not be aligned to any reference.
    Unmappable(String),
    IO(std::io::Error),
}

impl std::error::Error for RibosomeError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            RibosomeError::IO(e) => Some(e),
            _ => None,
        }
    }
}
impl GetCode for RibosomeError {
    fn get_code(&self) -> i32 {
        match self {
            RibosomeError::IO(e) => e.get_code(),
            RibosomeError::BlankFirstLine(_) => 66, // EX_NOINPUT
            _ => 1,
        }
    }
}

impl std::fmt::Display for RibosomeError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            RibosomeError::InvalidFastaFormat => {
                write!(f, "Invalid FASTA format. Header needs to be either ID or ID|ctype.")
            }
            RibosomeError::InvalidTsvFormat => write!(f, "Invalid TSV format: expected 2 or 3 tab-separated columns"),
            RibosomeError::BlankFirstLine(p) => write!(f, "Blank or empty first line: {}", p.display()),
            RibosomeError::NoCtype(s) => write!(f, "No ctype found in header: {s}"),
            RibosomeError::UnimplementedCtype(c) => write!(f, "Unimpelmented ctype '{c}'"),
            RibosomeError::InvalidSequence(s) => write!(f, "Invalid sequence, see header: {s}"),
            RibosomeError::Unmappable(id) => write!(f, "Query '{id}' could not be aligned to any reference"),
            RibosomeError::IO(_) => write!(f, "An underlying parse or IO failure occurred"),
        }
    }
}

impl From<std::io::Error> for RibosomeError {
    fn from(e: std::io::Error) -> Self {
        RibosomeError::IO(e)
    }
}

impl From<ErrorWithContext> for RibosomeError {
    fn from(e: ErrorWithContext) -> Self {
        RibosomeError::IO(std::io::Error::other(e))
    }
}
