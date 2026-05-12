//! Row structs, parsers, and display implementations for the genome sequence
//! file.

use crate::{
    outputs::PrecomputedGenomeData,
    toml::Formatting,
    tsv::{HADOOP_NULL, Nullable},
};
use csv::{Reader, ReaderBuilder};
use serde::Deserialize;
use serde_derive::Deserialize;
use std::{fmt::Display, fs::File, io::Read, path::Path};
use zoe::{
    data::err::ResultWithErrorContext,
    prelude::{DataOwned, Nucleotides, NucleotidesView},
};

/// The data in a single row of the genome sequence file.
#[derive(Clone, Eq, PartialEq, Hash, Debug)]
pub struct GenSeqRow {
    pub query_id:      String,
    pub ctype:         String,
    pub reference_id:  String,
    pub genome_id:     Option<String>,
    pub genome_length: usize,
    pub has_insertion: bool,
    pub genome_seq:    Nucleotides,
    pub genome_aln:    Nucleotides,
}

impl<'de> Deserialize<'de> for GenSeqRow {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>, {
        let GenSeqRowRaw {
            query_id,
            ctype,
            reference_id,
            genome_id,
            genome_length,
            has_insertion,
            genome_seq,
            genome_aln,
        } = GenSeqRowRaw::deserialize(deserializer)?;

        let genome_id = Nullable::from(genome_id).into_option();
        let genome_seq = Nucleotides::from(genome_seq);
        let genome_aln = Nucleotides::from(genome_aln);

        Ok(GenSeqRow {
            query_id,
            ctype,
            reference_id,
            genome_id,
            genome_length,
            has_insertion,
            genome_seq,
            genome_aln,
        })
    }
}

/// A helper struct for deserializing [`GenSeqRow`].
#[derive(Deserialize)]
struct GenSeqRowRaw {
    query_id:      String,
    ctype:         String,
    reference_id:  String,
    genome_id:     String,
    genome_length: usize,
    has_insertion: bool,
    genome_seq:    String,
    genome_aln:    String,
}

/// The data in a single row of the genome sequence file, with all fields
/// borrowed.
///
/// This is useful for writing a [`GenSeqRow`] record without needing to
/// clone/allocate each part.
#[derive(Clone, Eq, PartialEq, Hash, Debug)]
pub struct GenSeqRowView<'a> {
    pub query_id:        &'a str,
    pub ctype:           &'a str,
    pub reference_id:    &'a str,
    pub genome_id:       Option<&'a str>,
    pub genome_length:   usize,
    pub has_insertion:   bool,
    pub genome_seq:      NucleotidesView<'a>,
    pub genome_aln:      NucleotidesView<'a>,
    pub genome_aln_rpad: usize,
}

impl<'a> GenSeqRowView<'a> {
    /// Creates a new [`GenSeqRowView`] by extracting the relevant fields from
    /// the [`PrecomputedGenomeData`].
    pub fn new(
        genome: &'a PrecomputedGenomeData, query_id: &'a str, ctype: &'a str, reference_id: &'a str,
        formatting: &'a Formatting,
    ) -> Self {
        let genome_aln_rpad = if formatting.right_pad_gen {
            genome.trailing_ref_unaligned
        } else {
            0
        };

        Self {
            query_id,
            ctype,
            reference_id,
            genome_id: genome.genome_id.as_ref().map(AsRef::as_ref),
            genome_length: genome.genome_length,
            has_insertion: genome.has_insertion,
            genome_seq: genome.genome_seq.as_view(),
            genome_aln: genome.genome_aln.as_view(),
            genome_aln_rpad,
        }
    }
}

impl Display for GenSeqRow {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            concat!(
                "{query_id}\t{ctype}\t{reference_id}\t{genome_id}",
                "\t{genome_length}\t{has_insertion}\t{genome_seq}\t{genome_aln}",
            ),
            query_id = self.query_id,
            ctype = self.ctype,
            reference_id = self.reference_id,
            genome_id = self.genome_id.as_ref().map(AsRef::as_ref).unwrap_or(HADOOP_NULL),
            genome_length = self.genome_length,
            has_insertion = self.has_insertion,
            genome_seq = self.genome_seq,
            genome_aln = self.genome_aln,
        )
    }
}

impl Display for GenSeqRowView<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            concat!(
                "{query_id}\t{ctype}\t{reference_id}\t{genome_id}",
                "\t{genome_length}\t{has_insertion}\t{genome_seq}",
                "\t{genome_aln}{empty:.<genome_aln_rpad$}"
            ),
            query_id = self.query_id,
            ctype = self.ctype,
            reference_id = self.reference_id,
            genome_id = self.genome_id.as_ref().map(AsRef::as_ref).unwrap_or(HADOOP_NULL),
            genome_length = self.genome_length,
            has_insertion = self.has_insertion,
            genome_seq = self.genome_seq,
            genome_aln = self.genome_aln,
            empty = "",
            genome_aln_rpad = self.genome_aln_rpad,
        )
    }
}

/// A parser for the genome sequence file output by DAIS-ribosome.
pub struct GenSeqFileParser<R: Read> {
    reader: Reader<R>,
}

impl GenSeqFileParser<File> {
    /// Opens a new [`GenSeqFileParser`] from a provided `path`.
    ///
    /// ## Errors
    ///
    /// Any IO errors while opening the file are propagated with context.
    pub fn from_file(path: impl AsRef<Path>) -> std::io::Result<Self> {
        Ok(Self::from_readable(
            File::open(&path).with_path_context("Failed to open genome seq file:", path)?,
        ))
    }
}

impl<R: Read> GenSeqFileParser<R> {
    /// Creates a new [`GenSeqFileParser`] from a provided `readable` type.
    pub fn from_readable(readable: R) -> Self {
        Self {
            reader: ReaderBuilder::new().has_headers(false).delimiter(b'\t').from_reader(readable),
        }
    }
}

impl<R: Read> Iterator for GenSeqFileParser<R> {
    type Item = Result<GenSeqRow, csv::Error>;

    fn next(&mut self) -> Option<Self::Item> {
        self.reader.deserialize().next()
    }
}
