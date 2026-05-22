//! Row structs, parsers, and display implementations for the product sequence
//! file.

use crate::{
    outputs::ComputedProduct,
    ranges::{CdsCoord, InclusiveDisplay, parse_coords_inclusive},
    toml::Formatting,
    tsv::{HADOOP_NULL, Nullable},
};
use csv::{Reader, ReaderBuilder};
use serde::{Deserialize, de::Error};
use serde_derive::Deserialize;
use std::{fmt::Display, fs::File, io::Read, ops::Range, path::Path};
use zoe::{
    data::err::{ResultWithErrorContext, WithErrorContext},
    prelude::{AminoAcids, AminoAcidsView, DataOwned, Len, Nucleotides, NucleotidesView},
};

/// The data in a single row of the product sequence file.
#[derive(Clone, Eq, PartialEq, Hash, Debug)]
pub struct SeqRow {
    /// The ID of the query.
    pub query_id:          String,
    /// The compound type of the query.
    pub ctype:             String,
    /// The ID for the reference group which was aligned against.
    pub reference_id:      String,
    /// The protein product name (e.g., `HA`, `HA-signal`).
    pub product_name:      String,
    /// The MD5 hash of the cleaned amino acid sequence (variant hash), or
    /// `None` if no amino acid data remained after filtering.
    ///
    /// See [`ComputedProduct::variant_hash`].
    pub variant_hash:      Option<String>,
    /// The unaligned amino acid sequence for the protein (with insertions but
    /// no deletions).
    ///
    /// See [`ComputedProduct::aa_seq`].
    pub aa_seq:            Option<AminoAcids>,
    /// The aligned amino acid sequence for the protein (with `-` for deletions
    /// but no insertions).
    ///
    /// See [`ComputedProduct::aa_aln`].
    pub aa_aln:            Option<AminoAcids>,
    /// The SHA1 hash of the cleaned coding sequence, or `None` if no DNA data
    /// remained after filtering.
    ///
    /// See [`ComputedProduct::cds_id`].
    pub cds_id:            Option<String>,
    /// Whether any insertion exists in this product.
    ///
    /// See [`ComputedProduct::has_insertion`].
    pub has_insertion:     bool,
    /// Whether any insertion or deletion causes a frameshift (i.e., the length
    /// is not divisible by 3).
    ///
    /// See [`ComputedProduct::has_shift_indel`].
    pub has_shift_indel:   bool,
    /// The unaligned coding sequence for the protein (with insertions but no
    /// deletions).
    ///
    /// See [`ComputedProduct::cds_seq`].
    pub cds_seq:           Option<Nucleotides>,
    /// The aligned coding sequence for the protein (with `-` for deletions but
    /// no insertions).
    ///
    /// See [`ComputedProduct::cds_aln`].
    pub cds_aln:           Option<Nucleotides>,
    /// The coordinates within the original query that were used to form the
    /// unaligned `cds_seq`.
    ///
    /// See [`ComputedProduct::query_coords`].
    pub query_coordinates: Vec<Range<usize>>,
    /// The coding sequence coordinates corresponding to the
    /// `query_coordinates`.
    ///
    /// See [`ComputedProduct::cds_coords`].
    pub cds_coordinates:   Vec<CdsCoord>,
}

impl<'de> Deserialize<'de> for SeqRow {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>, {
        let SeqRowRaw {
            query_id,
            ctype,
            reference_id,
            product_name,
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

        let query_coordinates = parse_coords_inclusive::<Range<usize>>(&query_coordinates)
            .collect::<Result<Vec<_>, _>>()
            .map_err(|e| {
                D::Error::custom(e.with_context(format!("Failed to parse query coordinates: {query_coordinates}")))
            })?;

        let cds_coordinates = parse_coords_inclusive::<CdsCoord>(&cds_coordinates)
            .collect::<Result<Vec<_>, _>>()
            .map_err(|e| {
                D::Error::custom(e.with_context(format!("Failed to parse coding sequence coordinates: {cds_coordinates}")))
            })?;

        Ok(SeqRow {
            query_id,
            ctype,
            reference_id,
            product_name,
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
    product_name:      String,
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
    /// The ID of the query.
    pub query_id:          &'a str,
    /// The compound type of the query.
    pub ctype:             &'a str,
    /// The ID for the reference group which was aligned against.
    pub reference_id:      &'a str,
    /// The protein product name (e.g., `HA`, `HA-signal`).
    pub product_name:      &'a str,
    /// The MD5 hash of the cleaned amino acid sequence (variant hash), or
    /// `None` if no amino acid data remained after filtering.
    pub variant_hash:      Option<&'a str>,
    /// The unaligned amino acid sequence for the protein (with insertions but
    /// no deletions).
    ///
    /// See [`ComputedProduct::aa_seq`].
    pub aa_seq:            AminoAcidsView<'a>,
    /// The aligned amino acid sequence for the protein (with `-` for deletions
    /// but no insertions).
    ///
    /// See [`ComputedProduct::aa_aln`].
    pub aa_aln:            AminoAcidsView<'a>,
    /// The amount of right padding to apply to `aa_aln`.
    pub aa_aln_rpad:       usize,
    /// The SHA1 hash of the cleaned coding sequence, or `None` if no DNA data
    /// remained after filtering.
    pub cds_id:            Option<&'a str>,
    /// Whether any insertion exists in this product.
    pub has_insertion:     bool,
    /// Whether any insertion or deletion causes a frameshift (i.e., the length
    /// is not divisible by 3).
    pub has_shift_indel:   bool,
    /// The unaligned coding sequence for the protein (with insertions but no
    /// deletions).
    ///
    /// See [`ComputedProduct::cds_seq`].
    pub cds_seq:           NucleotidesView<'a>,
    /// The aligned coding sequence for the protein (with `-` for deletions but
    /// no insertions).
    ///
    /// See [`ComputedProduct::cds_aln`].
    pub cds_aln:           NucleotidesView<'a>,
    /// The amount of right padding to apply to `cds_aln`.
    pub cds_aln_rpad:      usize,
    /// The coordinates within the original query that were used to form the
    /// unaligned `cds_seq`.
    ///
    /// See [`ComputedProduct::query_coords`].
    pub query_coordinates: &'a [Range<usize>],
    /// The coding sequence coordinates corresponding to the
    /// `query_coordinates`.
    ///
    /// See [`ComputedProduct::cds_coords`].
    pub cds_coordinates:   &'a [CdsCoord],
}

impl<'a> SeqRowView<'a> {
    /// Creates a new [`SeqRowView`] by extracting the relevant fields from the
    /// [`ComputedProduct`].
    #[allow(unused_variables)]
    pub fn new(
        product: &'a ComputedProduct<'a>, query_id: &'a str, ctype: &'a str, reference_id: &'a str,
        formatting: &'a Formatting,
    ) -> Self {
        let cds_aln_rpad = if formatting.right_pad_cds && !product.cds_aln.is_empty() {
            product.trailing_cds_unaligned
        } else {
            0
        };
        let aa_aln_rpad = if formatting.right_pad_aa && !product.aa_aln.is_empty() {
            // Floor divide since any leftover bases are already accounted for with a partial codon
            product.trailing_cds_unaligned / 3
        } else {
            0
        };

        Self {
            query_id,
            ctype,
            reference_id,
            product_name: product.name,
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
                "{query_id}\t{ctype}\t{reference_id}\t{product_name}\t{vh}",
                "\t{aa_seq}\t{aa_aln}",
                "\t{cds_id}\t{ins}\t{shift}",
                "\t{cds_seq}\t{cds_aln}",
                "\t{query_coordinates}\t{cds_coordinates}"
            ),
            query_id = self.query_id,
            ctype = self.ctype,
            reference_id = self.reference_id,
            product_name = self.product_name,
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
            query_coordinates = self.query_coordinates.display_inclusive(),
            cds_coordinates = self.cds_coordinates.display_inclusive(),
        )
    }
}

impl Display for SeqRowView<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            concat!(
                "{query_id}\t{ctype}\t{reference_id}\t{product_name}\t{vh}",
                "\t{aa_seq}\t{aa_aln}{empty:.<aa_aln_rpad$}",
                "\t{cds_id}\t{ins}\t{shift}",
                "\t{cds_seq}\t{cds_aln}{empty:.<cds_aln_rpad$}",
                "\t{query_coordinates}\t{cds_coordinates}"
            ),
            query_id = self.query_id,
            ctype = self.ctype,
            reference_id = self.reference_id,
            product_name = self.product_name,
            vh = self.variant_hash.unwrap_or(HADOOP_NULL),
            aa_seq = Nullable(self.aa_seq),
            aa_aln = Nullable(&self.aa_aln),
            aa_aln_rpad = self.aa_aln_rpad,
            cds_id = self.cds_id.unwrap_or(HADOOP_NULL),
            ins = self.has_insertion,
            shift = self.has_shift_indel,
            cds_seq = Nullable(self.cds_seq),
            cds_aln = Nullable(self.cds_aln),
            cds_aln_rpad = self.cds_aln_rpad,
            query_coordinates = self.query_coordinates.display_inclusive(),
            cds_coordinates = self.cds_coordinates.display_inclusive(),
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
            File::open(&path).with_path_context("Failed to open seq file:", path)?,
        ))
    }
}

impl<R: Read> SeqFileParser<R> {
    /// Creates a new [`SeqFileParser`] from a provided `readable` type.
    pub fn from_readable(readable: R) -> Self {
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
