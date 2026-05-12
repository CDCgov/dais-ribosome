//! Error types for DAIS-ribosome.

use zoe::data::err::{ErrorWithContext, GetCode};

/// An error holding a compound type that was not implemented for the given
/// module.
#[derive(Clone, Eq, PartialEq, Hash, Debug)]
pub struct UnimplementedCtype(pub String);

impl std::error::Error for UnimplementedCtype {}
impl GetCode for UnimplementedCtype {}

impl std::fmt::Display for UnimplementedCtype {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Unimpelmented ctype '{}'", self.0)
    }
}

impl From<String> for UnimplementedCtype {
    fn from(ctype: String) -> Self {
        Self(ctype)
    }
}

/// A minimal enum containing error variants which must be matched on in
/// DAIS-ribosome. All other errors should use [`RibosomeError::Io`].
#[derive(Debug)]
pub enum RibosomeError {
    UnimplementedCtype(UnimplementedCtype),
    Io(std::io::Error),
}

impl std::error::Error for RibosomeError {
    #[inline]
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            RibosomeError::Io(e) => e.source(),
            _ => None,
        }
    }
}

impl GetCode for RibosomeError {
    #[inline]
    fn get_code(&self) -> i32 {
        match self {
            RibosomeError::UnimplementedCtype(e) => e.get_code(),
            RibosomeError::Io(e) => e.get_code(),
        }
    }
}

impl std::fmt::Display for RibosomeError {
    #[inline]
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            RibosomeError::UnimplementedCtype(e) => write!(f, "{e}"),
            RibosomeError::Io(e) => write!(f, "{e}"),
        }
    }
}

impl From<&str> for RibosomeError {
    #[inline]
    fn from(value: &str) -> Self {
        RibosomeError::Io(std::io::Error::other(value))
    }
}

impl From<String> for RibosomeError {
    #[inline]
    fn from(value: String) -> Self {
        RibosomeError::Io(std::io::Error::other(value))
    }
}

impl From<std::io::Error> for RibosomeError {
    #[inline]
    fn from(e: std::io::Error) -> Self {
        RibosomeError::Io(e)
    }
}

impl From<ErrorWithContext> for RibosomeError {
    #[inline]
    fn from(e: ErrorWithContext) -> Self {
        RibosomeError::Io(std::io::Error::other(e))
    }
}
