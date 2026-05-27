//! Row structs, parsers, and display implementations for the product sequence
//! file.

use crate::{
    outputs::{ComputedProduct, MaybeComputedProduct},
    ranges::{CdsCoord, InclusiveDisplay, parse_coords_inclusive},
    toml::Formatting,
    tsv::Nullable,
};
use csv::{Reader, ReaderBuilder};
use serde::{Deserialize, de::Error};
use std::{fmt::Display, fs::File, io::Read, ops::Range, path::Path};
use zoe::{
    data::err::{ResultWithErrorContext, WithErrorContext},
    prelude::{AminoAcids, AminoAcidsView, DataOwned, Len, Nucleotides, NucleotidesView},
};

pub enum SeqRow {
    Data(SeqData),
    Empty(EmptySeqRow),
    Deleted(DeletedSeqRow),
}

/// The data in a single row of the product sequence file.
#[derive(Clone, Eq, PartialEq, Hash, Debug)]
pub struct EmptySeqRow {
    /// The ID of the query.
    pub query_id:     String,
    /// The compound type of the query.
    pub ctype:        String,
    /// The ID for the reference group which was aligned against.
    pub reference_id: String,
    /// The protein product name (e.g., `HA`, `HA-signal`).
    pub product_name: String,
}

/// The data in a single row of the product sequence file.
#[derive(Clone, Eq, PartialEq, Hash, Debug)]
pub struct DeletedSeqRow {
    /// The ID of the query.
    pub query_id:     String,
    /// The compound type of the query.
    pub ctype:        String,
    /// The ID for the reference group which was aligned against.
    pub reference_id: String,
    /// The protein product name (e.g., `HA`, `HA-signal`).
    pub product_name: String,
    /// The aligned amino acid sequence for the protein (with `-` for deletions
    /// but no insertions).
    ///
    /// See [`ComputedProduct::aa_aln`].
    pub aa_aln:       AminoAcids,
    /// The aligned coding sequence for the protein (with `-` for deletions but
    /// no insertions).
    ///
    /// See [`ComputedProduct::cds_aln`].
    pub cds_aln:      Nucleotides,
}

/// The data in a single row of the product sequence file.
#[derive(Clone, Eq, PartialEq, Hash, Debug)]
pub struct SeqData {
    /// The ID of the query.
    pub query_id:          String,
    /// The compound type of the query.
    pub ctype:             String,
    /// The ID for the reference group which was aligned against.
    pub reference_id:      String,
    /// The protein product name (e.g., `HA`, `HA-signal`).
    pub product_name:      String,
    /// The MD5 hash of the cleaned amino acid sequence (variant hash).
    ///
    /// See [`ComputedProduct::variant_hash`].
    pub variant_hash:      String,
    /// The unaligned amino acid sequence for the protein (with insertions but
    /// no deletions).
    ///
    /// See [`ComputedProduct::aa_seq`].
    pub aa_seq:            AminoAcids,
    /// The aligned amino acid sequence for the protein (with `-` for deletions
    /// but no insertions).
    ///
    /// See [`ComputedProduct::aa_aln`].
    pub aa_aln:            AminoAcids,
    /// The SHA1 hash of the cleaned coding sequence, or `None` if no DNA data
    /// remained after filtering.
    ///
    /// See [`ComputedProduct::cds_id`].
    pub cds_id:            String,
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
    pub cds_seq:           Nucleotides,
    /// The aligned coding sequence for the protein (with `-` for deletions but
    /// no insertions).
    ///
    /// See [`ComputedProduct::cds_aln`], except that this may also contain
    /// right padding.
    pub cds_aln:           Nucleotides,
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
    /// Deserializes a [`SeqRow`] into one of the three variants.
    ///
    /// ## Errors
    ///
    /// The combination of null fields must be valid (conforming to one of the
    /// three variants). The coordinates must successfully be parsed, otherwise
    /// an error with context is yielded.
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

        let Some(aa_aln) = aa_aln.into_option() else {
            let missing_aa_aln_requirements = [
                ("aa_seq", aa_seq, "aa_aln"),
                ("variant_hash", variant_hash, "aa_seq"),
                ("cds_aln", cds_aln, "aa_aln"),
                ("cds_seq", cds_seq, "cds_aln"),
                ("cds_id", cds_id, "cds_seq"),
                ("query_coordinates", query_coordinates, "cds_seq"),
                ("cds_coordinates", cds_coordinates, "cds_seq"),
            ];

            for (field_name, field, requires_name) in missing_aa_aln_requirements {
                if let Some(field) = field.into_option() {
                    return Err(D::Error::custom(format!(
                        "{field_name} requires {requires_name} to be present. Found null for {requires_name} and the following for {field_name}: {field}"
                    )));
                }
            }

            return Ok(Self::Empty(EmptySeqRow {
                query_id,
                ctype,
                reference_id,
                product_name,
            }));
        };

        let aa_aln = AminoAcids::from(aa_aln);

        let Some(cds_aln) = cds_aln.into_option() else {
            return Err(D::Error::custom(format!(
                "aa_aln requires cds_aln to be present. Found null for cds_aln and the following for aa_aln: {aa_aln}"
            )));
        };

        let cds_aln = Nucleotides::from(cds_aln);

        let Some(aa_seq) = aa_seq.into_option() else {
            let missing_aa_seq_requirements = [
                ("variant_hash", variant_hash, "aa_seq"),
                ("cds_seq", cds_seq, "aa_seq"),
                ("cds_id", cds_id, "cds_seq"),
                ("query_coordinates", query_coordinates, "cds_seq"),
                ("cds_coordinates", cds_coordinates, "cds_seq"),
            ];

            for (field_name, field, requires_name) in missing_aa_seq_requirements {
                if let Some(field) = field.into_option() {
                    return Err(D::Error::custom(format!(
                        "{field_name} requires {requires_name} to be present. Found null for {requires_name} and the following for {field_name}: {field}"
                    )));
                }
            }

            return Ok(Self::Deleted(DeletedSeqRow {
                query_id,
                ctype,
                reference_id,
                product_name,
                aa_aln,
                cds_aln,
            }));
        };

        let aa_seq = AminoAcids::from(aa_seq);

        let Some(cds_seq) = cds_seq.into_option() else {
            return Err(D::Error::custom(format!(
                "aa_seq requires cds_seq to be present. Found null for cds_seq and the following for aa_seq: {aa_seq}"
            )));
        };

        let cds_seq = Nucleotides::from(cds_seq);

        let Some(variant_hash) = variant_hash.into_option() else {
            return Err(D::Error::custom(format!(
                "aa_seq requires variant_hash to be present. Found null for variant_hash and the following for aa_seq: {aa_seq}"
            )));
        };

        let Some(cds_id) = cds_id.into_option() else {
            return Err(D::Error::custom(format!(
                "cds_seq requires cds_id to be present. Found null for cds_id and the following for cds_seq: {cds_seq}"
            )));
        };

        let Some(query_coordinates) = query_coordinates.into_option() else {
            return Err(D::Error::custom(format!(
                "cds_seq requires query_coordinates to be present. Found null for query_coordinates and the following for cds_seq: {cds_seq}"
            )));
        };

        let query_coordinates = parse_coords_inclusive::<Range<usize>>(&query_coordinates)
            .collect::<Result<Vec<_>, _>>()
            .map_err(|e| {
                D::Error::custom(e.with_context(format!("Failed to parse query coordinates: {query_coordinates}")))
            })?;

        let Some(cds_coordinates) = cds_coordinates.into_option() else {
            return Err(D::Error::custom(format!(
                "cds_seq requires cds_coordinates to be present. Found null for cds_coordinates and the following for cds_seq: {cds_seq}"
            )));
        };

        let cds_coordinates = parse_coords_inclusive::<CdsCoord>(&cds_coordinates)
            .collect::<Result<Vec<_>, _>>()
            .map_err(|e| {
                D::Error::custom(e.with_context(format!("Failed to parse coding sequence coordinates: {cds_coordinates}")))
            })?;

        Ok(SeqRow::Data(SeqData {
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
        }))
    }
}

/// A helper struct for deserializing [`SeqRow`].
#[derive(Deserialize)]
struct SeqRowRaw {
    query_id:          String,
    ctype:             String,
    reference_id:      String,
    product_name:      String,
    variant_hash:      Nullable<String>,
    aa_seq:            Nullable<String>,
    aa_aln:            Nullable<String>,
    cds_id:            Nullable<String>,
    has_insertion:     bool,
    has_shift_indel:   bool,
    cds_seq:           Nullable<String>,
    cds_aln:           Nullable<String>,
    query_coordinates: Nullable<String>,
    cds_coordinates:   Nullable<String>,
}

pub enum SeqRowView<'a> {
    Data(SeqDataView<'a>),
    Empty(EmptySeqRowView<'a>),
    Deleted(DeletedSeqRowView<'a>),
}

impl<'a> SeqRowView<'a> {
    /// Creates a new [`SeqRowView`] by extracting the relevant fields from the
    /// [`MaybeComputedProduct`].
    pub fn new(
        product: &'a MaybeComputedProduct<'a>, query_id: &'a str, ctype: &'a str, reference_id: &'a str,
        formatting: &'a Formatting,
    ) -> Self {
        match product {
            MaybeComputedProduct::Ok(product) => {
                SeqRowView::Data(SeqDataView::new(product, query_id, ctype, reference_id, formatting))
            }
            MaybeComputedProduct::Empty(product) => SeqRowView::Empty(EmptySeqRowView {
                query_id,
                ctype,
                reference_id,
                product_name: product.name,
            }),
            MaybeComputedProduct::Deleted(product) => {
                let cds_aln_rpad = if formatting.right_pad_cds && !product.cds_aln.is_empty() {
                    product.cds_aln_rpad
                } else {
                    0
                };
                let aa_aln_rpad = if formatting.right_pad_aa && !product.aa_aln.is_empty() {
                    product.aa_aln_rpad()
                } else {
                    0
                };

                SeqRowView::Deleted(DeletedSeqRowView {
                    query_id,
                    ctype,
                    reference_id,
                    product_name: product.name,
                    aa_aln: product.aa_aln.as_view(),
                    aa_aln_rpad,
                    cds_aln: product.cds_aln.as_view(),
                    cds_aln_rpad,
                })
            }
        }
    }
}

/// The data in a single row of the product sequence file.
#[derive(Clone, Eq, PartialEq, Hash, Debug)]
pub struct EmptySeqRowView<'a> {
    /// The ID of the query.
    pub query_id:     &'a str,
    /// The compound type of the query.
    pub ctype:        &'a str,
    /// The ID for the reference group which was aligned against.
    pub reference_id: &'a str,
    /// The protein product name (e.g., `HA`, `HA-signal`).
    pub product_name: &'a str,
}

/// The data in a single row of the product sequence file.
#[derive(Clone, Eq, PartialEq, Hash, Debug)]
pub struct DeletedSeqRowView<'a> {
    /// The ID of the query.
    pub query_id:     &'a str,
    /// The compound type of the query.
    pub ctype:        &'a str,
    /// The ID for the reference group which was aligned against.
    pub reference_id: &'a str,
    /// The protein product name (e.g., `HA`, `HA-signal`).
    pub product_name: &'a str,
    /// The aligned amino acid sequence for the protein (with `-` for deletions
    /// but no insertions).
    ///
    /// See [`ComputedProduct::aa_aln`].
    pub aa_aln:       AminoAcidsView<'a>,
    /// The amount of right padding to apply to `aa_aln`.
    pub aa_aln_rpad:  usize,
    /// The aligned coding sequence for the protein (with `-` for deletions but
    /// no insertions).
    ///
    /// See [`ComputedProduct::cds_aln`].
    pub cds_aln:      NucleotidesView<'a>,
    /// The amount of right padding to apply to `cds_aln`.
    pub cds_aln_rpad: usize,
}

/// The data in a single row of the product sequence file, with all fields
/// borrowed.
///
/// This is useful for writing a [`SeqRow`] record without needing to
/// clone/allocate each part.
#[derive(Clone, Eq, PartialEq, Hash, Debug)]
pub struct SeqDataView<'a> {
    /// The ID of the query.
    pub query_id:          &'a str,
    /// The compound type of the query.
    pub ctype:             &'a str,
    /// The ID for the reference group which was aligned against.
    pub reference_id:      &'a str,
    /// The protein product name (e.g., `HA`, `HA-signal`).
    pub product_name:      &'a str,
    /// The MD5 hash of the cleaned amino acid sequence (variant hash).
    pub variant_hash:      &'a str,
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
    /// The amount of right padding to apply to `aa_aln` when displaying it.
    ///
    /// For example, if [`Formatting::right_pad_aa`] is false or the product is
    /// empty, then this is set to 0 by [`SeqRowView::new`].
    pub aa_aln_rpad:       usize,
    /// The SHA1 hash of the cleaned coding sequence, or `None` if no DNA data
    /// remained after filtering.
    pub cds_id:            &'a str,
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
    /// The amount of right padding to apply to `cds_aln` when displaying it.
    ///
    /// For example, if [`Formatting::right_pad_cds`] is false or the product is
    /// empty, then this is set to 0 by [`SeqRowView::new`].
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

impl<'a> SeqDataView<'a> {
    /// Creates a new [`SeqDataView`] by extracting the relevant fields from the
    /// [`ComputedProduct`].
    pub fn new(
        product: &'a ComputedProduct<'a>, query_id: &'a str, ctype: &'a str, reference_id: &'a str,
        formatting: &'a Formatting,
    ) -> Self {
        let cds_aln_rpad = if formatting.right_pad_cds && !product.cds_aln.is_empty() {
            product.cds_aln_rpad
        } else {
            0
        };
        let aa_aln_rpad = if formatting.right_pad_aa && !product.aa_aln.is_empty() {
            product.aa_aln_rpad()
        } else {
            0
        };

        Self {
            query_id,
            ctype,
            reference_id,
            product_name: product.name,
            variant_hash: &product.variant_hash,
            aa_seq: product.aa_seq.as_view(),
            aa_aln: product.aa_aln.as_view(),
            aa_aln_rpad,
            cds_id: &product.cds_id,
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
        match self {
            SeqRow::Data(row) => row.fmt(f),
            SeqRow::Empty(row) => row.fmt(f),
            SeqRow::Deleted(row) => row.fmt(f),
        }
    }
}

impl Display for SeqData {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let SeqData {
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
        } = self;

        let display = SeqRowDisplay {
            query_id,
            ctype,
            reference_id,
            product_name,
            variant_hash: Nullable(variant_hash),
            aa_seq: Nullable(aa_seq.as_view()),
            aa_aln: Nullable(aa_aln.as_view()),
            aa_aln_rpad: 0,
            cds_id: Nullable(cds_id),
            has_insertion: *has_insertion,
            has_shift_indel: *has_shift_indel,
            cds_seq: Nullable(cds_seq.as_view()),
            cds_aln: Nullable(cds_aln.as_view()),
            cds_aln_rpad: 0,
            query_coordinates: Nullable(query_coordinates),
            cds_coordinates: Nullable(cds_coordinates),
        };

        write!(f, "{display}")
    }
}

impl Display for EmptySeqRow {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let EmptySeqRow {
            query_id,
            ctype,
            reference_id,
            product_name,
        } = self;

        let display = SeqRowDisplay {
            query_id,
            ctype,
            reference_id,
            product_name,
            ..Default::default()
        };

        write!(f, "{display}")
    }
}

impl Display for DeletedSeqRow {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let DeletedSeqRow {
            query_id,
            ctype,
            reference_id,
            product_name,
            aa_aln,
            cds_aln,
        } = self;

        let display = SeqRowDisplay {
            query_id,
            ctype,
            reference_id,
            product_name,
            aa_aln: Nullable(aa_aln.as_view()),
            cds_aln: Nullable(cds_aln.as_view()),
            ..Default::default()
        };

        write!(f, "{display}")
    }
}

impl Display for SeqRowView<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SeqRowView::Data(row) => row.fmt(f),
            SeqRowView::Empty(row) => row.fmt(f),
            SeqRowView::Deleted(row) => row.fmt(f),
        }
    }
}

impl Display for SeqDataView<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let SeqDataView {
            query_id,
            ctype,
            reference_id,
            product_name,
            variant_hash,
            aa_seq,
            aa_aln,
            aa_aln_rpad,
            cds_id,
            has_insertion,
            has_shift_indel,
            cds_seq,
            cds_aln,
            cds_aln_rpad,
            query_coordinates,
            cds_coordinates,
        } = *self;

        let display = SeqRowDisplay {
            query_id,
            ctype,
            reference_id,
            product_name,
            variant_hash: Nullable(variant_hash),
            aa_seq: Nullable(aa_seq),
            aa_aln: Nullable(aa_aln),
            aa_aln_rpad,
            cds_id: Nullable(cds_id),
            has_insertion,
            has_shift_indel,
            cds_seq: Nullable(cds_seq),
            cds_aln: Nullable(cds_aln),
            cds_aln_rpad,
            query_coordinates: Nullable(query_coordinates),
            cds_coordinates: Nullable(cds_coordinates),
        };

        write!(f, "{display}")
    }
}

impl Display for EmptySeqRowView<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let EmptySeqRowView {
            query_id,
            ctype,
            reference_id,
            product_name,
        } = self;

        let display = SeqRowDisplay {
            query_id,
            ctype,
            reference_id,
            product_name,
            ..Default::default()
        };

        write!(f, "{display}")
    }
}

impl Display for DeletedSeqRowView<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let DeletedSeqRowView {
            query_id,
            ctype,
            reference_id,
            product_name,
            aa_aln,
            aa_aln_rpad,
            cds_aln,
            cds_aln_rpad,
        } = *self;

        let display = SeqRowDisplay {
            query_id,
            ctype,
            reference_id,
            product_name,
            aa_aln: Nullable(aa_aln),
            aa_aln_rpad,
            cds_aln: Nullable(cds_aln),
            cds_aln_rpad,
            ..Default::default()
        };

        write!(f, "{display}")
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

/// A helper struct for displaying sequence-related TSV output.
///
/// All nullable fields are represented as `Option`. This struct is used to
/// remove redundant [`Display`] implementations and ensure greater correctness.
#[derive(Default)]
struct SeqRowDisplay<'a> {
    query_id:          &'a str,
    ctype:             &'a str,
    reference_id:      &'a str,
    product_name:      &'a str,
    variant_hash:      Nullable<&'a str>,
    aa_seq:            Nullable<AminoAcidsView<'a>>,
    aa_aln:            Nullable<AminoAcidsView<'a>>,
    aa_aln_rpad:       usize,
    cds_id:            Nullable<&'a str>,
    has_insertion:     bool,
    has_shift_indel:   bool,
    cds_seq:           Nullable<NucleotidesView<'a>>,
    cds_aln:           Nullable<NucleotidesView<'a>>,
    cds_aln_rpad:      usize,
    query_coordinates: Nullable<&'a [Range<usize>]>,
    cds_coordinates:   Nullable<&'a [CdsCoord]>,
}

impl Display for SeqRowDisplay<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            concat!(
                "{query_id}\t{ctype}\t{reference_id}\t{product_name}",
                "\t{variant_hash}\t{aa_seq}\t{aa_aln}{empty:.<aa_aln_rpad$}",
                "\t{cds_id}\t{has_insertion}\t{has_shift_indel}",
                "\t{cds_seq}\t{cds_aln}{empty:.<cds_aln_rpad$}",
                "\t{query_coordinates}\t{cds_coordinates}"
            ),
            query_id = self.query_id,
            ctype = self.ctype,
            reference_id = self.reference_id,
            product_name = self.product_name,
            variant_hash = self.variant_hash,
            aa_seq = self.aa_seq,
            aa_aln = self.aa_aln,
            aa_aln_rpad = self.aa_aln_rpad,
            cds_id = self.cds_id,
            has_insertion = self.has_insertion,
            has_shift_indel = self.has_shift_indel,
            cds_seq = self.cds_seq,
            cds_aln = self.cds_aln,
            cds_aln_rpad = self.cds_aln_rpad,
            query_coordinates = self
                .query_coordinates
                .display_null_or_else(|coords| coords.display_inclusive()),
            cds_coordinates = self.cds_coordinates.display_null_or_else(|coords| coords.display_inclusive()),
            empty = ""
        )
    }
}
