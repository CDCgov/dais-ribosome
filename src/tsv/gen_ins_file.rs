//! Row structs, parsers, and display implementations for the genome insertion
//! file.

use crate::outputs::ComputedGenomeInsertion;
use csv::{Reader, ReaderBuilder};
use serde::Deserialize;
use serde_derive::Deserialize;
use std::{fmt::Display, fs::File, io::Read, path::Path};
use zoe::{
    data::err::ResultWithErrorContext,
    prelude::{DataOwned, Nucleotides, NucleotidesView},
};

/// The data in a single row of the genome insertion file.
#[derive(Clone, Eq, PartialEq, Hash, Debug)]
pub struct GenInsRow {
    pub query_id:     String,
    pub ctype:        String,
    pub reference_id: String,
    pub nt_pos:       usize,
    pub inserted_nts: Nucleotides,
}

impl<'de> Deserialize<'de> for GenInsRow {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>, {
        let GenInsRowRaw {
            query_id,
            ctype,
            reference_id,
            nt_pos,
            inserted_nts,
        } = GenInsRowRaw::deserialize(deserializer)?;

        let inserted_nts = Nucleotides::from(inserted_nts);

        Ok(GenInsRow {
            query_id,
            ctype,
            reference_id,
            nt_pos,
            inserted_nts,
        })
    }
}

/// A helper struct for deserializing [`GenInsRow`].
#[derive(Deserialize)]
struct GenInsRowRaw {
    query_id:     String,
    ctype:        String,
    reference_id: String,
    nt_pos:       usize,
    inserted_nts: String,
}

/// The data in a single row of the genome insertion file, with all fields
/// borrowed.
///
/// This is useful for writing a [`GenInsRow`] record without needing to
/// clone/allocate each part.
#[derive(Clone, Eq, PartialEq, Hash, Debug)]
pub struct GenInsRowView<'a> {
    pub query_id:     &'a str,
    pub ctype:        &'a str,
    pub reference_id: &'a str,
    pub nt_pos:       usize,
    pub inserted_nts: NucleotidesView<'a>,
}

impl<'a> GenInsRowView<'a> {
    /// Creates a new [`GenInsRowView`] by extracting the relevant fields from
    /// the [`ComputedGenomeInsertion`].
    pub fn new(insertion: &'a ComputedGenomeInsertion, query_id: &'a str, ctype: &'a str, reference_id: &'a str) -> Self {
        Self {
            query_id,
            ctype,
            reference_id,
            nt_pos: insertion.upstream_nt_pos,
            inserted_nts: insertion.inserted_nucleotides.as_view(),
        }
    }
}

impl Display for GenInsRow {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}\t{}\t{}\t{}\t{}",
            self.query_id, self.ctype, self.reference_id, self.nt_pos, self.inserted_nts,
        )
    }
}

impl Display for GenInsRowView<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}\t{}\t{}\t{}\t{}",
            self.query_id, self.ctype, self.reference_id, self.nt_pos, self.inserted_nts,
        )
    }
}

/// A parser for the genome insertion file output by DAIS-ribosome.
pub struct GenInsFileParser<R: Read> {
    reader: Reader<R>,
}

impl GenInsFileParser<File> {
    /// Opens a new [`GenInsFileParser`] from a provided `path`.
    ///
    /// ## Errors
    ///
    /// Any IO errors while opening the file are propagated with context.
    pub fn from_file(path: impl AsRef<Path>) -> std::io::Result<Self> {
        Ok(Self::from_readable(
            File::open(&path).with_path_context("Failed to open genome ins file:", path)?,
        ))
    }
}

impl<R: Read> GenInsFileParser<R> {
    /// Creates a new [`GenInsFileParser`] from a provided `readable` type.
    pub fn from_readable(readable: R) -> Self {
        Self {
            reader: ReaderBuilder::new().has_headers(false).delimiter(b'\t').from_reader(readable),
        }
    }
}

impl<R: Read> Iterator for GenInsFileParser<R> {
    type Item = Result<GenInsRow, csv::Error>;

    fn next(&mut self) -> Option<Self::Item> {
        self.reader.deserialize().next()
    }
}
