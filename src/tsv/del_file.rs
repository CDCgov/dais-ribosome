//! Row structs, parsers, and display implementations for the product deletion
//! file.

use crate::{
    data::products::{ComputedDeletion, ComputedProduct},
    tsv::{HADOOP_NULL, Nullable},
};
use csv::{Reader, ReaderBuilder};
use serde::Deserialize;
use serde_derive::Deserialize;
use std::{fmt::Display, fs::File, io::Read, path::Path};
use zoe::data::err::ResultWithErrorContext;

/// The data in a single row of the product deletion file.
#[derive(Clone, Eq, PartialEq, Hash, Debug)]
pub struct DelRow {
    pub query_id:      String,
    pub ctype:         String,
    pub reference_id:  String,
    pub protein:       String,
    pub variant_hash:  Option<String>,
    pub del_aa_start:  usize,
    pub del_aa_end:    usize,
    pub del_aa_len:    usize,
    pub in_frame:      bool,
    pub cds_id:        Option<String>,
    pub del_cds_start: usize,
    pub del_cds_end:   usize,
    pub del_cds_len:   usize,
}

impl<'de> Deserialize<'de> for DelRow {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>, {
        let DelRowRaw {
            query_id,
            ctype,
            reference_id,
            protein,
            variant_hash,
            del_aa_start,
            del_aa_end,
            del_aa_len,
            in_frame,
            cds_id,
            del_cds_start,
            del_cds_end,
            del_cds_len,
        } = DelRowRaw::deserialize(deserializer)?;

        let variant_hash = Nullable::from(variant_hash).into_option();
        let cds_id = Nullable::from(cds_id).into_option();

        Ok(DelRow {
            query_id,
            ctype,
            reference_id,
            protein,
            variant_hash,
            del_aa_start,
            del_aa_end,
            del_aa_len,
            in_frame,
            cds_id,
            del_cds_start,
            del_cds_end,
            del_cds_len,
        })
    }
}

/// A helper struct for deserializing [`DelRow`].
#[derive(Deserialize)]
struct DelRowRaw {
    query_id:      String,
    ctype:         String,
    reference_id:  String,
    protein:       String,
    variant_hash:  String,
    del_aa_start:  usize,
    del_aa_end:    usize,
    del_aa_len:    usize,
    in_frame:      bool,
    cds_id:        String,
    del_cds_start: usize,
    del_cds_end:   usize,
    del_cds_len:   usize,
}

/// The data in a single row of the product deletion file, with all fields
/// borrowed.
///
/// This is useful for writing a [`DelRow`] record without needing to
/// clone/allocate each part.
#[derive(Clone, Eq, PartialEq, Hash, Debug)]
pub struct DelRowView<'a> {
    pub query_id:      &'a str,
    pub ctype:         &'a str,
    pub reference_id:  &'a str,
    pub protein:       &'a str,
    pub variant_hash:  Option<&'a str>,
    pub del_aa_start:  usize,
    pub del_aa_end:    usize,
    pub del_aa_len:    usize,
    pub in_frame:      bool,
    pub cds_id:        Option<&'a str>,
    pub del_cds_start: usize,
    pub del_cds_end:   usize,
    pub del_cds_len:   usize,
}

impl<'a> DelRowView<'a> {
    /// Creates a new [`DelRowView`] by extracting the relevant fields from the
    /// [`ComputedDeletion`] and [`ComputedProduct`].
    pub fn new(
        deletion: &'a ComputedDeletion, product: &'a ComputedProduct<'a>, query_id: &'a str, ctype: &'a str,
        reference_id: &'a str,
    ) -> Self {
        Self {
            query_id,
            ctype,
            reference_id,
            protein: product.product_name,
            variant_hash: product.variant_hash.as_ref().map(AsRef::as_ref),
            del_aa_start: deletion.del_aa_start,
            del_aa_end: deletion.del_aa_end,
            del_aa_len: deletion.del_aa_len,
            in_frame: deletion.in_frame,
            cds_id: product.cds_id.as_ref().map(AsRef::as_ref),
            del_cds_start: deletion.del_cds_start,
            del_cds_end: deletion.del_cds_end,
            del_cds_len: deletion.del_cds_len,
        }
    }
}

impl Display for DelRow {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.query_id,
            self.ctype,
            self.reference_id,
            self.protein,
            self.variant_hash.as_ref().map(AsRef::as_ref).unwrap_or(HADOOP_NULL),
            self.del_aa_start,
            self.del_aa_end,
            self.del_aa_len,
            self.in_frame,
            self.cds_id.as_ref().map(AsRef::as_ref).unwrap_or(HADOOP_NULL),
            self.del_cds_start,
            self.del_cds_end,
            self.del_cds_len,
        )
    }
}

impl Display for DelRowView<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.query_id,
            self.ctype,
            self.reference_id,
            self.protein,
            self.variant_hash.as_ref().map(AsRef::as_ref).unwrap_or(HADOOP_NULL),
            self.del_aa_start,
            self.del_aa_end,
            self.del_aa_len,
            self.in_frame,
            self.cds_id.as_ref().map(AsRef::as_ref).unwrap_or(HADOOP_NULL),
            self.del_cds_start,
            self.del_cds_end,
            self.del_cds_len,
        )
    }
}

/// A parser for the product deletion file output by DAIS-ribosome.
pub struct DelFileParser<R: Read> {
    reader: Reader<R>,
}

impl DelFileParser<File> {
    /// Opens a new [`DelFileParser`] from a provided `path`.
    ///
    /// ## Errors
    ///
    /// Any IO errors while opening the file are propagated with context.
    pub fn from_file(path: impl AsRef<Path>) -> std::io::Result<Self> {
        Ok(Self::from_readable(
            File::open(&path).with_file_context("Failed to open del file:", path)?,
        ))
    }
}

impl<R: Read> DelFileParser<R> {
    /// Creates a new [`DelFileParser`] from a provided `readable` type.
    pub fn from_readable(readable: R) -> Self {
        Self {
            reader: ReaderBuilder::new().has_headers(false).delimiter(b'\t').from_reader(readable),
        }
    }
}

impl<R: Read> Iterator for DelFileParser<R> {
    type Item = Result<DelRow, csv::Error>;

    fn next(&mut self) -> Option<Self::Item> {
        self.reader.deserialize().next()
    }
}
