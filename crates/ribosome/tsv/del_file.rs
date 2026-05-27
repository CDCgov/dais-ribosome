//! Row structs, parsers, and display implementations for the product deletion
//! file.

use crate::{
    outputs::{ComputedDeletion, ComputedProduct, DeletedProduct},
    tsv::{HADOOP_NULL, Nullable},
};
use csv::{Reader, ReaderBuilder};
use serde::Deserialize;
use std::{fmt::Display, fs::File, io::Read, path::Path};
use zoe::data::err::ResultWithErrorContext;

/// The data in a single row of the product deletion file.
#[derive(Clone, Eq, PartialEq, Hash, Debug)]
pub struct DelRow {
    /// The ID of the query.
    pub query_id:      String,
    /// The compound type of the query.
    pub ctype:         String,
    /// The ID for the reference group which was aligned against.
    pub reference_id:  String,
    /// The protein product name (e.g., `HA`, `HA-signal`).
    pub product_name:  String,
    /// The variant hash of the product.
    ///
    /// See [`ComputedProduct::variant_hash`].
    pub variant_hash:  Option<String>,
    /// The start position of the deletion in amino acid coordinates (1-based,
    /// inclusive).
    ///
    /// See [`ComputedDeletion::del_aa_start`].
    pub del_aa_start:  usize,
    /// The end position of the deletion in amino acid coordinates (1-based,
    /// inclusive).
    ///
    /// See [`ComputedDeletion::del_aa_end`].
    pub del_aa_end:    usize,
    /// The deletion length in amino acids.
    ///
    /// See [`ComputedDeletion::del_aa_len`].
    pub del_aa_len:    usize,
    /// Whether deletion is in-frame (both the CDS start position and length are
    /// multiples of 3).
    ///
    /// See [`ComputedDeletion::in_frame`].
    pub in_frame:      bool,
    /// The ID of the coding sequence.
    ///
    /// See [`ComputedProduct::cds_id`].
    pub cds_id:        Option<String>,
    /// The start position of the deletion in coding sequence coordinates
    /// (1-based, inclusive).
    ///
    /// See [`ComputedDeletion::del_cds_start`].
    pub del_cds_start: usize,
    /// The end position of the deletion in coding sequence coordinates
    /// (1-based, inclusive).
    ///
    /// See [`ComputedDeletion::del_cds_end`].
    pub del_cds_end:   usize,
    /// The deletion length in nucleotides.
    ///
    /// See [`ComputedDeletion::del_cds_len`].
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
            product_name,
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

        Ok(DelRow {
            query_id,
            ctype,
            reference_id,
            product_name,
            variant_hash: variant_hash.into_option(),
            del_aa_start,
            del_aa_end,
            del_aa_len,
            in_frame,
            cds_id: cds_id.into_option(),
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
    product_name:  String,
    variant_hash:  Nullable<String>,
    del_aa_start:  usize,
    del_aa_end:    usize,
    del_aa_len:    usize,
    in_frame:      bool,
    cds_id:        Nullable<String>,
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
    /// The ID of the query.
    pub query_id:      &'a str,
    /// The compound type of the query.
    pub ctype:         &'a str,
    /// The ID for the reference group which was aligned against.
    pub reference_id:  &'a str,
    /// The protein product name (e.g., `HA`, `HA-signal`).
    pub product_name:  &'a str,
    /// The variant hash of the product.
    ///
    /// See [`ComputedProduct::variant_hash`].
    pub variant_hash:  Option<&'a str>,
    /// The start position of the deletion in amino acid coordinates (1-based,
    /// inclusive).
    ///
    /// See [`ComputedDeletion::del_aa_start`].
    pub del_aa_start:  usize,
    /// The end position of the deletion in amino acid coordinates (1-based,
    /// inclusive).
    ///
    /// See [`ComputedDeletion::del_aa_end`].
    pub del_aa_end:    usize,
    /// The deletion length in amino acids.
    ///
    /// See [`ComputedDeletion::del_aa_len`].
    pub del_aa_len:    usize,
    /// Whether deletion is in-frame (both the CDS start position and length are
    /// multiples of 3).
    ///
    /// See [`ComputedDeletion::in_frame`].
    pub in_frame:      bool,
    /// The ID of the coding sequence.
    ///
    /// See [`ComputedProduct::cds_id`].
    pub cds_id:        Option<&'a str>,
    /// The start position of the deletion in coding sequence coordinates
    /// (1-based, inclusive).
    ///
    /// See [`ComputedDeletion::del_cds_start`].
    pub del_cds_start: usize,
    /// The end position of the deletion in coding sequence coordinates
    /// (1-based, inclusive).
    ///
    /// See [`ComputedDeletion::del_cds_end`].
    pub del_cds_end:   usize,
    /// The deletion length in nucleotides.
    ///
    /// See [`ComputedDeletion::del_cds_len`].
    pub del_cds_len:   usize,
}

impl<'a> DelRowView<'a> {
    /// Creates a new [`DelRowView`] by extracting the relevant fields from the
    /// [`ComputedDeletion`] and [`ComputedProduct`].
    ///
    /// The `variant_hash` and `cds_id` fields will be `Some`.
    pub fn new(
        deletion: &'a ComputedDeletion, product: &'a ComputedProduct<'a>, query_id: &'a str, ctype: &'a str,
        reference_id: &'a str,
    ) -> Self {
        Self {
            query_id,
            ctype,
            reference_id,
            product_name: product.name,
            variant_hash: Some(&product.variant_hash),
            del_aa_start: deletion.del_aa_start,
            del_aa_end: deletion.del_aa_end,
            del_aa_len: deletion.del_aa_len,
            in_frame: deletion.in_frame,
            cds_id: Some(&product.cds_id),
            del_cds_start: deletion.del_cds_start,
            del_cds_end: deletion.del_cds_end,
            del_cds_len: deletion.del_cds_len,
        }
    }

    /// Creates a new [`DelRowView`] by extracting the relevant fields from the
    /// [`DeletedProduct`].
    ///
    /// The `variant_hash` and `cds_id` fields will be `None`.
    pub fn from_deleted_product(
        product: &'a DeletedProduct<'a>, query_id: &'a str, ctype: &'a str, reference_id: &'a str,
    ) -> Self {
        let deletion = &product.deletion;

        Self {
            query_id,
            ctype,
            reference_id,
            product_name: product.name,
            variant_hash: None,
            del_aa_start: deletion.del_aa_start,
            del_aa_end: deletion.del_aa_end,
            del_aa_len: deletion.del_aa_len,
            in_frame: deletion.in_frame,
            cds_id: None,
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
            self.product_name,
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
            self.product_name,
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
            File::open(&path).with_path_context("Failed to open del file:", path)?,
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
