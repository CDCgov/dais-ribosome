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
    pub query_id:     String,
    pub ctype:        String,
    pub reference_id: String,
    pub protein:      String,
    pub aa_pos:       usize,
    pub inserted_nts: Nucleotides,
    pub inserted_aas: AminoAcids,
    pub nt_pos:       usize,
    pub codon_shift:  usize,
}

impl<'de> Deserialize<'de> for InsRow {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>, {
        let InsRowRaw {
            query_id,
            ctype,
            reference_id,
            protein,
            aa_pos,
            inserted_nts,
            inserted_aas,
            nt_pos,
            codon_shift,
        } = InsRowRaw::deserialize(deserializer)?;

        let inserted_nts = Nucleotides::from(inserted_nts);
        let inserted_aas = AminoAcids::from(inserted_aas);

        Ok(InsRow {
            query_id,
            ctype,
            reference_id,
            protein,
            aa_pos,
            inserted_nts,
            inserted_aas,
            nt_pos,
            codon_shift,
        })
    }
}

/// A helper struct for deserializing [`InsRow`].
#[derive(Deserialize)]
struct InsRowRaw {
    query_id:     String,
    ctype:        String,
    reference_id: String,
    protein:      String,
    aa_pos:       usize,
    inserted_nts: String,
    inserted_aas: String,
    nt_pos:       usize,
    codon_shift:  usize,
}

/// The data in a single row of the product insertion file, with all fields
/// borrowed.
///
/// This is useful for writing a [`InsRow`] record without needing to
/// clone/allocate each part.
#[derive(Clone, Eq, PartialEq, Hash, Debug)]
pub struct InsRowView<'a> {
    pub query_id:     &'a str,
    pub ctype:        &'a str,
    pub reference_id: &'a str,
    pub protein:      &'a str,
    pub aa_pos:       usize,
    pub inserted_nts: NucleotidesView<'a>,
    pub inserted_aas: AminoAcidsView<'a>,
    pub nt_pos:       usize,
    pub codon_shift:  usize,
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
            protein: product.product_name,
            aa_pos: insertion.upstream_aa_pos,
            inserted_nts: insertion.inserted_nucleotides.as_view(),
            inserted_aas: insertion.inserted_residues.as_view(),
            nt_pos: insertion.upstream_nt_pos,
            codon_shift: insertion.codon_shift,
        }
    }
}

impl Display for InsRowView<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        #[cfg(feature = "regression-testing")]
        let nts = self.inserted_nts.to_string().to_lowercase();

        #[cfg(not(feature = "regression-testing"))]
        let nts = &self.inserted_nts;

        write!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.query_id,
            self.ctype,
            self.reference_id,
            self.protein,
            self.aa_pos,
            nts,
            self.inserted_aas,
            self.nt_pos,
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
