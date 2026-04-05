//! Row structs, parsers, and display implementations for the product sequence
//! file.

use crate::{
    config::Formatting,
    data::products::ComputedProduct,
    tsv::{HADOOP_NULL, Nullable},
};
use csv::{Reader, ReaderBuilder};
use serde::Deserialize;
use serde_derive::Deserialize;
use std::{fmt::Display, fs::File, io::Read, path::Path};
use zoe::{
    data::err::ResultWithErrorContext,
    prelude::{AminoAcids, AminoAcidsView, DataOwned, Len, Nucleotides, NucleotidesView},
};

/// The data in a single row of the product sequence file.
#[derive(Clone, Eq, PartialEq, Hash, Debug)]
pub struct SeqRow {
    pub query_id:          String,
    pub ctype:             String,
    pub reference_id:      String,
    pub protein:           String,
    pub variant_hash:      Option<String>,
    pub aa_seq:            Option<AminoAcids>,
    pub aa_aln:            Option<AminoAcids>,
    pub cds_id:            Option<String>,
    pub has_insertion:     bool,
    pub has_shift_indel:   bool,
    pub cds_seq:           Option<Nucleotides>,
    pub cds_aln:           Option<Nucleotides>,
    pub query_coordinates: String,
    pub cds_coordinates:   String,
}

impl<'de> Deserialize<'de> for SeqRow {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>, {
        let SeqRowRaw {
            query_id,
            ctype,
            reference_id,
            protein,
            variant_hash,
            aa_seq,
            aa_aln,
            cds_id,
            has_insertion,
            has_shift_indel,
            cds_seq,
            cds_aln,
            query_coordinates,
            cds_coordinates,
        } = SeqRowRaw::deserialize(deserializer)?;

        let variant_hash = Nullable::from(variant_hash).into_option();
        let aa_seq = Nullable::from(aa_seq).into_option();
        let aa_aln = Nullable::from(aa_aln).into_option();
        let cds_id = Nullable::from(cds_id).into_option();
        let cds_seq = Nullable::from(cds_seq).into_option();
        let cds_aln = Nullable::from(cds_aln).into_option();

        Ok(SeqRow {
            query_id,
            ctype,
            reference_id,
            protein,
            variant_hash,
            aa_seq,
            aa_aln,
            cds_id,
            has_insertion,
            has_shift_indel,
            cds_seq,
            cds_aln,
            query_coordinates,
            cds_coordinates,
        })
    }
}

/// A helper struct for deserializing [`SeqRow`].
#[derive(Deserialize)]
struct SeqRowRaw {
    query_id:          String,
    ctype:             String,
    reference_id:      String,
    protein:           String,
    variant_hash:      String,
    aa_seq:            String,
    aa_aln:            String,
    cds_id:            String,
    has_insertion:     bool,
    has_shift_indel:   bool,
    cds_seq:           String,
    cds_aln:           String,
    query_coordinates: String,
    cds_coordinates:   String,
}

/// The data in a single row of the product sequence file, with all fields
/// borrowed.
///
/// This is useful for writing a [`SeqRow`] record without needing to
/// clone/allocate each part.
#[derive(Clone, Eq, PartialEq, Hash, Debug)]
pub struct SeqRowView<'a> {
    pub query_id:          &'a str,
    pub ctype:             &'a str,
    pub reference_id:      &'a str,
    pub protein:           &'a str,
    pub variant_hash:      Option<&'a str>,
    pub aa_seq:            AminoAcidsView<'a>,
    pub aa_aln:            AminoAcidsView<'a>,
    pub aa_aln_rpad:       usize,
    pub cds_id:            Option<&'a str>,
    pub has_insertion:     bool,
    pub has_shift_indel:   bool,
    pub cds_seq:           NucleotidesView<'a>,
    pub cds_aln:           NucleotidesView<'a>,
    pub cds_aln_rpad:      usize,
    pub query_coordinates: &'a str,
    pub cds_coordinates:   &'a str,
}

impl<'a> SeqRowView<'a> {
    /// Creates a new [`SeqRowView`] by extracting the relevant fields from the
    /// [`ComputedProduct`].
    #[allow(unused_variables)]
    pub fn new(
        product: &'a ComputedProduct<'a>, query_id: &'a str, ctype: &'a str, reference_id: &'a str,
        formatting: &'a Formatting,
    ) -> Self {
        // Regression: always pad CDS, never pad AA, synthesize empty AA alignments.
        #[cfg(feature = "regression-testing")]
        let cds_aln_rpad = product.trailing_cds_unaligned;
        #[cfg(feature = "regression-testing")]
        let aa_aln_rpad = if product.aa_aln.is_empty() {
            (product.cds_aln.len() + product.trailing_cds_unaligned) / 3
        } else {
            0
        };

        // Normal: skip padding when there is no data
        #[cfg(not(feature = "regression-testing"))]
        let cds_aln_rpad = if formatting.right_pad_cds && !product.cds_aln.is_empty() {
            product.trailing_cds_unaligned
        } else {
            0
        };
        #[cfg(not(feature = "regression-testing"))]
        let aa_aln_rpad = if formatting.right_pad_aa && !product.aa_aln.is_empty() {
            product.trailing_cds_unaligned / 3
        } else {
            0
        };

        Self {
            query_id,
            ctype,
            reference_id,
            protein: product.product_name,
            variant_hash: product.variant_hash.as_ref().map(AsRef::as_ref),
            aa_seq: product.aa_seq.as_view(),
            aa_aln: product.aa_aln.as_view(),
            aa_aln_rpad,
            cds_id: product.cds_id.as_ref().map(AsRef::as_ref),
            has_insertion: product.has_insertion,
            has_shift_indel: product.has_shift_indel,
            cds_seq: product.cds_seq.as_view(),
            cds_aln: product.cds_aln.as_view(),
            cds_aln_rpad,
            query_coordinates: product.query_coords.as_ref(),
            cds_coordinates: product.cds_coords.as_ref(),
        }
    }
}

impl Display for SeqRow {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            concat!(
                "{query_id}\t{ctype}\t{reference_id}\t{protein}\t{vh}",
                "\t{aa_seq}\t{aa_aln}",
                "\t{cds_id}\t{ins}\t{shift}",
                "\t{cds_seq}\t{cds_aln}",
                "\t{query_coordinates}\t{cds_coordinates}"
            ),
            query_id = self.query_id,
            ctype = self.ctype,
            reference_id = self.reference_id,
            protein = self.protein,
            vh = self.variant_hash.as_ref().map(AsRef::as_ref).unwrap_or(HADOOP_NULL),
            aa_seq = self
                .aa_seq
                .as_ref()
                .map_or(AminoAcidsView::from(HADOOP_NULL), |s| s.as_view()),
            aa_aln = self
                .aa_aln
                .as_ref()
                .map_or(AminoAcidsView::from(HADOOP_NULL), |s| s.as_view()),
            cds_id = self.cds_id.as_ref().map(AsRef::as_ref).unwrap_or(HADOOP_NULL),
            ins = self.has_insertion,
            shift = self.has_shift_indel,
            cds_seq = self
                .cds_seq
                .as_ref()
                .map_or(NucleotidesView::from(HADOOP_NULL), |s| s.as_view()),
            cds_aln = self
                .cds_aln
                .as_ref()
                .map_or(NucleotidesView::from(HADOOP_NULL), |s| s.as_view()),
            query_coordinates = Nullable(&self.query_coordinates),
            cds_coordinates = Nullable(&self.cds_coordinates),
        )
    }
}

impl Display for SeqRowView<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        #[cfg(feature = "regression-testing")]
        let aa_aln = self.aa_aln;
        #[cfg(not(feature = "regression-testing"))]
        let aa_aln = Nullable(&computed_product.aa_aln);

        write!(
            f,
            concat!(
                "{query_id}\t{ctype}\t{reference_id}\t{protein}\t{vh}",
                "\t{aa_seq}\t{aa_aln}{empty:.<aa_aln_rpad$}",
                "\t{cds_id}\t{ins}\t{shift}",
                "\t{cds_seq}\t{cds_aln}{empty:.<cds_aln_rpad$}",
                "\t{query_coordinates}\t{cds_coordinates}"
            ),
            query_id = self.query_id,
            ctype = self.ctype,
            reference_id = self.reference_id,
            protein = self.protein,
            vh = self.variant_hash.unwrap_or(HADOOP_NULL),
            aa_seq = Nullable(self.aa_seq),
            aa_aln = aa_aln,
            aa_aln_rpad = self.aa_aln_rpad,
            cds_id = self.cds_id.unwrap_or(HADOOP_NULL),
            ins = self.has_insertion,
            shift = self.has_shift_indel,
            cds_seq = Nullable(self.cds_seq),
            cds_aln = self.cds_aln,
            cds_aln_rpad = self.cds_aln_rpad,
            query_coordinates = Nullable(self.query_coordinates),
            cds_coordinates = Nullable(self.cds_coordinates),
            empty = ""
        )
    }
}

/// A parser for the product sequence file output by DAIS-ribosome.
pub struct SeqFileParser<R: Read> {
    reader: Reader<R>,
}

impl SeqFileParser<File> {
    /// Opens a new [`SeqFileParser`] from a provided `path`.
    ///
    /// ## Errors
    ///
    /// Any IO errors while opening the file are propagated with context.
    pub fn from_file(path: impl AsRef<Path>) -> std::io::Result<Self> {
        Ok(Self::from_readable(
            File::open(&path).with_file_context("Failed to open seq file:", path)?,
        ))
    }
}

impl<R: Read> SeqFileParser<R> {
    /// Creates a new [`SeqFileParser`] from a provided `readable` type.
    pub fn from_readable(readable: R) -> Self {
        // TODO: Should we automatically detect headers? What process would add
        // headers? What would they be? MIRA-oxide had a case with headers.
        Self {
            reader: ReaderBuilder::new().has_headers(false).delimiter(b'\t').from_reader(readable),
        }
    }
}

impl<R: Read> Iterator for SeqFileParser<R> {
    type Item = Result<SeqRow, csv::Error>;

    fn next(&mut self) -> Option<Self::Item> {
        self.reader.deserialize().next()
    }
}
