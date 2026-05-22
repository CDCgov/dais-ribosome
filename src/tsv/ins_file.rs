//! Row structs, parsers, and display implementations for the product insertion
//! file.

use crate::outputs::{ComputedInsertion, ComputedProduct};
use csv::{Reader, ReaderBuilder};
use serde::Deserialize;
use serde_derive::Deserialize;
use std::{fmt::Display, fs::File, io::Read, path::Path};
use zoe::{
    data::err::ResultWithErrorContext,
    prelude::{AminoAcids, AminoAcidsView, DataOwned, Nucleotides, NucleotidesView},
};

/// The data in a single row of the product insertion file.
#[derive(Clone, Eq, PartialEq, Hash, Debug)]
pub struct InsRow {
    /// The ID of the query.
    pub query_id:        String,
    /// The compound type of the query.
    pub ctype:           String,
    /// The ID for the reference group which was aligned against.
    pub reference_id:    String,
    /// The protein product name (e.g., `HA`, `HA-signal`).
    pub product_name:    String,
    /// The upstream amino acid position (1-based), which is the position
    /// _after_ which the insertion occurs.
    ///
    /// See [`ComputedInsertion::upstream_aa_pos`].
    pub upstream_aa_pos: usize,
    /// The inserted nucleotides.
    ///
    /// See [`ComputedInsertion::inserted_nt`].
    pub inserted_nt:     Nucleotides,
    /// A direct translation of `inserted_nt` to amino acids.
    ///
    /// See [`ComputedInsertion::inserted_aa`].
    pub inserted_aa:     AminoAcids,
    /// The upstream nucleotide position (1-based), which is the position
    /// _after_ which the insertion occurs.
    ///
    /// See [`ComputedInsertion::upstream_nt_pos`].
    pub upstream_nt_pos: usize,
    /// The codon shift of the insertion.
    ///
    /// See [`ComputedInsertion::codon_shift`].
    pub codon_shift:     usize,
}

impl<'de> Deserialize<'de> for InsRow {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>, {
        let InsRowRaw {
            query_id,
            ctype,
            reference_id,
            product_name,
            upstream_aa_pos,
            inserted_nt,
            inserted_aa,
            upstream_nt_pos,
            codon_shift,
        } = InsRowRaw::deserialize(deserializer)?;

        let inserted_nt = Nucleotides::from(inserted_nt);
        let inserted_aa = AminoAcids::from(inserted_aa);

        Ok(InsRow {
            query_id,
            ctype,
            reference_id,
            product_name,
            upstream_aa_pos,
            inserted_nt,
            inserted_aa,
            upstream_nt_pos,
            codon_shift,
        })
    }
}

/// A helper struct for deserializing [`InsRow`].
#[derive(Deserialize)]
struct InsRowRaw {
    query_id:        String,
    ctype:           String,
    reference_id:    String,
    product_name:    String,
    upstream_aa_pos: usize,
    inserted_nt:     String,
    inserted_aa:     String,
    upstream_nt_pos: usize,
    codon_shift:     usize,
}

/// The data in a single row of the product insertion file, with all fields
/// borrowed.
///
/// This is useful for writing a [`InsRow`] record without needing to
/// clone/allocate each part.
#[derive(Clone, Eq, PartialEq, Hash, Debug)]
pub struct InsRowView<'a> {
    /// The ID of the query.
    pub query_id:        &'a str,
    /// The compound type of the query.
    pub ctype:           &'a str,
    /// The ID for the reference group which was aligned against.
    pub reference_id:    &'a str,
    /// The protein product name (e.g., `HA`, `HA-signal`).
    pub product_name:    &'a str,
    /// The upstream amino acid position (1-based), which is the position
    /// _after_ which the insertion occurs.
    ///
    /// See [`ComputedInsertion::upstream_aa_pos`].
    pub upstream_aa_pos: usize,
    /// The inserted nucleotides.
    ///
    /// See [`ComputedInsertion::inserted_nt`].
    pub inserted_nt:     NucleotidesView<'a>,
    /// A direct translation of `inserted_nt` to amino acids.
    ///
    /// See [`ComputedInsertion::inserted_aa`].
    pub inserted_aa:     AminoAcidsView<'a>,
    /// The upstream nucleotide position (1-based), which is the position
    /// _after_ which the insertion occurs.
    ///
    /// See [`ComputedInsertion::upstream_nt_pos`].
    pub upstream_nt_pos: usize,
    /// The codon shift of the insertion.
    ///
    /// See [`ComputedInsertion::codon_shift`].
    pub codon_shift:     usize,
}

impl<'a> InsRowView<'a> {
    /// Creates a new [`InsRowView`] by extracting the relevant fields from the
    /// [`ComputedInsertion`] and [`ComputedProduct`].
    pub fn new(
        insertion: &'a ComputedInsertion, product: &'a ComputedProduct, query_id: &'a str, ctype: &'a str,
        reference_id: &'a str,
    ) -> Self {
        Self {
            query_id,
            ctype,
            reference_id,
            product_name: product.name,
            upstream_aa_pos: insertion.upstream_aa_pos,
            inserted_nt: insertion.inserted_nt.as_view(),
            inserted_aa: insertion.inserted_aa.as_view(),
            upstream_nt_pos: insertion.upstream_nt_pos,
            codon_shift: insertion.codon_shift,
        }
    }
}

impl Display for InsRowView<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.query_id,
            self.ctype,
            self.reference_id,
            self.product_name,
            self.upstream_aa_pos,
            self.inserted_nt,
            self.inserted_aa,
            self.upstream_nt_pos,
            self.codon_shift,
        )
    }
}

/// A parser for the product insertion file output by DAIS-ribosome.
pub struct InsFileParser<R: Read> {
    reader: Reader<R>,
}

impl InsFileParser<File> {
    /// Opens a new [`InsFileParser`] from a provided `path`.
    ///
    /// ## Errors
    ///
    /// Any IO errors while opening the file are propagated with context.
    pub fn from_file(path: impl AsRef<Path>) -> std::io::Result<Self> {
        Ok(Self::from_readable(
            File::open(&path).with_path_context("Failed to open ins file:", path)?,
        ))
    }
}

impl<R: Read> InsFileParser<R> {
    /// Creates a new [`InsFileParser`] from a provided `readable` type.
    pub fn from_readable(readable: R) -> Self {
        Self {
            reader: ReaderBuilder::new().has_headers(false).delimiter(b'\t').from_reader(readable),
        }
    }
}

impl<R: Read> Iterator for InsFileParser<R> {
    type Item = Result<InsRow, csv::Error>;

    fn next(&mut self) -> Option<Self::Item> {
        self.reader.deserialize().next()
    }
}
