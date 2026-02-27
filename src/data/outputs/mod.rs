mod cds;
mod genome;
mod output;

pub(crate) use cds::*;
pub(crate) use genome::*;
pub use output::*;
use std::fmt::{Display, Formatter, Result};
use zoe::prelude::{AminoAcids, Nucleotides};

use crate::data::products::incremental_products::Coords;

#[derive(Copy, Clone)]
pub(crate) struct Nullable<T>(pub(crate) T);

impl<T> Nullable<T> {
    pub fn map<U, F>(self, f: F) -> Nullable<U>
    where
        F: FnOnce(T) -> U, {
        Nullable(f(self.0))
    }
}

impl<T: Display + AsRef<[u8]>> Nullable<T> {
    fn fmt_if_nonempty(&self, f: &mut Formatter<'_>) -> Result {
        if self.0.as_ref().is_empty() {
            f.write_str(HADOOP_NULL)
        } else {
            write!(f, "{}", self.0)
        }
    }
}

pub(crate) const HADOOP_NULL: &str = "\\N";

impl Display for Nullable<&Coords> {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        self.map(|c| c.as_str()).fmt_if_nonempty(f)
    }
}

impl Display for Nullable<&Nucleotides> {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        self.fmt_if_nonempty(f)
    }
}

impl Display for Nullable<&AminoAcids> {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        self.fmt_if_nonempty(f)
    }
}

impl Display for Nullable<&Option<String>> {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        match &self.0 {
            Some(val) => write!(f, "{val}"),
            None => f.write_str(HADOOP_NULL),
        }
    }
}
