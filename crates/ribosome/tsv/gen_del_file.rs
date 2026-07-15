//! Row structs, parsers, and display implementations for the genome deletion
//! file.

use crate::data::ranges::DeletionRange;
use csv::{Reader, ReaderBuilder};
use serde::Deserialize;
use std::{fmt::Display, fs::File, io::Read, path::Path};
use zoe::data::err::ResultWithErrorContext;

/// The data in a single row of the genome deletion file.
#[derive(Clone, Eq, PartialEq, Hash, Debug, Deserialize)]
pub struct GenDelRow {
    /// The ID of the query.
    pub query_id:     String,
    /// The compound type of the query.
    pub ctype:        String,
    /// The ID for the reference group which was aligned against.
    pub reference_id: String,
    /// The start position of the deletion in the nucleotide sequence (1-based,
    /// inclusive).
    pub del_start:    usize,
    /// The end position of the deletion in the nucleotide sequence (1-based,
    /// inclusive).
    pub del_end:      usize,
    /// The deletion length in nucleotides.
    pub del_len:      usize,
}

/// The data in a single row of the genome deletion file, with all fields
/// borrowed.
///
/// This is useful for writing a [`GenDelRow`] record without needing to
/// clone/allocate each part.
#[derive(Clone, Eq, PartialEq, Hash, Debug)]
pub struct GenDelRowView<'a> {
    /// The ID of the query.
    pub query_id:     &'a str,
    /// The compound type of the query.
    pub ctype:        &'a str,
    /// The ID for the reference group which was aligned against.
    pub reference_id: &'a str,
    /// The start position of the deletion in the nucleotide sequence (1-based,
    /// inclusive).
    pub del_start:    usize,
    /// The end position of the deletion in the nucleotide sequence (1-based,
    /// inclusive).
    pub del_end:      usize,
    /// The deletion length in nucleotides.
    pub del_len:      usize,
}

impl<'a> GenDelRowView<'a> {
    /// Creates a new [`GenDelRowView`] by extracting the relevant fields from
    /// the [`DeletionRange`].
    ///
    /// ## Validity
    ///
    /// `query_id`, `ctype`, and `reference_id` cannot contain tabs.
    pub fn new(deletion: &'a DeletionRange, query_id: &'a str, ctype: &'a str, reference_id: &'a str) -> GenDelRowView<'a> {
        Self {
            query_id,
            ctype,
            reference_id,
            del_start: deletion.ref_range.start + 1,
            del_end: deletion.ref_range.end,
            del_len: deletion.ref_range.len(),
        }
    }
}

impl Display for GenDelRow {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{}",
            self.query_id, self.ctype, self.reference_id, self.del_start, self.del_end, self.del_len,
        )
    }
}

impl Display for GenDelRowView<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{}",
            self.query_id, self.ctype, self.reference_id, self.del_start, self.del_end, self.del_len,
        )
    }
}

/// A parser for the genome deletion file output by DAIS-ribosome.
pub struct GenDelFileParser<R: Read> {
    reader: Reader<R>,
}

impl GenDelFileParser<File> {
    /// Opens a new [`GenDelFileParser`] from a provided `path`.
    ///
    /// ## Errors
    ///
    /// Any IO errors while opening the file are propagated with context.
    pub fn from_file(path: impl AsRef<Path>) -> std::io::Result<Self> {
        Ok(Self::from_readable(
            File::open(&path).with_path_context("Failed to open genome del file:", path)?,
        ))
    }
}

impl<R: Read> GenDelFileParser<R> {
    /// Creates a new [`GenDelFileParser`] from a provided `readable` type.
    pub fn from_readable(readable: R) -> Self {
        Self {
            reader: ReaderBuilder::new().has_headers(false).delimiter(b'\t').from_reader(readable),
        }
    }
}

impl<R: Read> Iterator for GenDelFileParser<R> {
    type Item = Result<GenDelRow, csv::Error>;

    fn next(&mut self) -> Option<Self::Item> {
        self.reader.deserialize().next()
    }
}
