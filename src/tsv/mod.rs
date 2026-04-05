//! Structs for formatting and parsing DAIS-ribosome output in TSV format, as
//! well as functions to write the entire TSV files from ribosome [`outputs`].
//!
//! For writing TSV output, the [`Writers`] type is provided which can be used
//! with [`write_product_output`] and [`write_genome_output`]. These functions
//! take [`RibosomeOutput`] and automatically handle materializing the outputs
//! and formatting the TSV records.
//!
//! For parsing TSV files, a dedicated parser (implementing [`Iterator`]) is
//! provided for each file type.
//!
//! For more fine-grained control, each of the six files has an owned struct
//! (useful for parsing) and a view struct (useful for writing without
//! re-allocating).
//!
//! [`outputs`]: crate::outputs
//! [`RibosomeOutput`]: crate::outputs::RibosomeOutput

mod del_file;
mod gen_del_file;
mod gen_ins_file;
mod gen_seq_file;
mod ins_file;
mod seq_file;
mod writing;

pub use del_file::*;
pub use gen_del_file::*;
pub use gen_ins_file::*;
pub use gen_seq_file::*;
pub use ins_file::*;
pub use seq_file::*;
pub use writing::*;

use std::fmt::Display;

#[derive(Copy, Clone)]
pub(crate) struct Nullable<T>(pub(crate) T);

impl<T: Display + AsRef<[u8]>> Nullable<T> {
    pub fn into_option(self) -> Option<T> {
        (!self.0.as_ref().is_empty()).then_some(self.0)
    }
}

pub(crate) const HADOOP_NULL: &str = "\\N";

impl<T: AsRef<[u8]> + Display> Display for Nullable<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.0.as_ref().is_empty() {
            f.write_str(HADOOP_NULL)
        } else {
            write!(f, "{}", self.0)
        }
    }
}

impl<T: From<String>> From<String> for Nullable<T> {
    fn from(mut value: String) -> Self {
        if value == HADOOP_NULL {
            value = String::new()
        };

        Nullable(T::from(value))
    }
}
