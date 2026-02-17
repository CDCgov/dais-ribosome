mod cds;
mod genome;
mod output;

pub(crate) use cds::*;
pub(crate) use genome::*;
pub use output::*;
use std::fmt::{Display, Formatter, Result};

pub(crate) struct Nullable<T>(pub(crate) T);
pub(crate) const HADOOP_NULL: &str = "\\N";

impl<T: Display + AsRef<[u8]>> Display for Nullable<T> {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        if self.0.as_ref().is_empty() {
            f.write_str(HADOOP_NULL)
        } else {
            write!(f, "{}", self.0)
        }
    }
}
