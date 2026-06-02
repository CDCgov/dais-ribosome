//! Multiplexed query input: auto-detects FASTA vs. TSV by peeking at the
//! first byte, then yields [`QueryInfo`]s through a unified iterator.
//! The file is opened only once.

use crate::{ClassificationStrategy, log::time_stamp};
use dais_ribosome::{NoNucleotides, QueryRecord};
use sswsort::{ClassificationResult, Strand};
use std::{
    error::Error,
    fmt::Display,
    fs::File,
    io::{BufRead, BufReader, Lines},
    path::Path,
};
use zoe::{
    data::{
        err::{GetCode, ResultWithErrorContext},
        fasta::FastaSeq,
    },
    define_whichever,
    prelude::*,
    unwrap_or_return_some_err,
};

/// Extracts the ctype from a [`ClassificationResult`], performing any logging
/// or filtering based on `rules`.
///
/// If the strand of the classification is [`Strand::Minus`], then the sequence
/// is converted to reverse complement as well.
fn handle_classification(
    classification: &ClassificationResult, id: &str, sequence: &mut Vec<u8>, verbose: bool,
) -> Option<String> {
    let (taxon, strand) = match classification {
        ClassificationResult::Unrecognizable { best_score } => {
            if verbose {
                if let Some(best_score) = best_score {
                    time_stamp(
                        &format!("The sequence for the following ID was unrecognizable with a score of {best_score}: {id}"),
                        true,
                    );
                } else {
                    time_stamp(&format!("The sequence for the following ID was unrecognizable: {id}"), true);
                }
            }
            return None;
        }
        ClassificationResult::Classification {
            taxon,
            best_score: _,
            strand,
        } => (taxon, strand),
        ClassificationResult::Chimeric { taxa } => {
            if verbose {
                let mut warning =
                    format!("The sequence for the following ID was detected as chimeric: {id}. Possible taxa: ");

                for part in taxa.iter().cloned().intersperse(", ") {
                    warning.push_str(part);
                }

                time_stamp(&warning, true);
            }

            return None;
        }
        ClassificationResult::UnusuallyLong {
            taxon,
            best_score: _,
            strand: _,
        } => {
            if verbose {
                time_stamp(
                    &format!("The sequence for the following ID was detected as unusually long for taxa {taxon}: {id}"),
                    true,
                );
            }

            return None;
        }
        ClassificationResult::Unresolvable { taxa, best_score } => {
            if verbose {
                let mut warning = format!(
                    "The sequence for the following ID was detected as unresolvable with a score of {best_score}: {id}. Possible taxa: "
                );

                for part in taxa.iter().cloned().intersperse(", ") {
                    warning.push_str(part);
                }

                time_stamp(&warning, true);
            }

            return None;
        }
    };

    match strand {
        Strand::Plus => {}
        Strand::Minus => {
            NucleotidesViewMut::from(sequence).make_reverse_complement();
        }
    }

    Some(taxon.to_string())
}

/// The information from reading a FASTA or TSV input, possibly with a compound
/// type.
pub struct QueryInfo {
    /// The ID of the query
    id:       String,
    /// The unsanitized query sequence
    sequence: Vec<u8>,
    /// The compound type, if provided
    ctype:    Option<String>,
}

impl QueryInfo {
    /// Converts the [`QueryInfo`] into a [`QueryRecord`] by classifying the
    /// ctype using the given [`ClassificationStrategy`] if needed. If the
    /// reverse strand is aligned against, then the reverse complement of the
    /// sequence is taken.
    ///
    /// ## Errors
    ///
    /// If `classification` is `None` and the `ctype` field of the query is
    /// missing, then [`NoCtype`] is returned.
    pub fn classify_and_prepare(
        mut self, classification: &Option<ClassificationStrategy>, verbose: bool,
    ) -> Result<Option<QueryRecord>, NoCtype> {
        let ctype = match (self.ctype, classification) {
            (Some(ctype), _) => ctype,
            (None, Some(ClassificationStrategy::SswSort(module))) => {
                let classification = module.classify(&self.sequence);
                let Some(ctype) = handle_classification(&classification, &self.id, &mut self.sequence, verbose) else {
                    return Ok(None);
                };
                ctype
            }
            (None, Some(ClassificationStrategy::Default(default_ctype))) => default_ctype.clone(),
            (None, None) => return Err(NoCtype { id: self.id }),
        };

        match QueryRecord::new(self.id, self.sequence, ctype) {
            Ok(record) => Ok(Some(record)),
            Err(NoNucleotides { id }) => {
                if verbose {
                    time_stamp(&format!("A sequence contained no unaligned DNA data. See ID: {id}"), true);
                }
                Ok(None)
            }
        }
    }
}

/// An error caused by the absense of a ctype within an input file.
#[derive(Clone, Eq, PartialEq, Hash, Debug)]
pub struct NoCtype {
    /// The ID of the record in the file.
    pub id: String,
}

impl Display for NoCtype {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "A query was missing a ctype: {id}", id = self.id)
    }
}

impl Error for NoCtype {}
impl GetCode for NoCtype {}

/// A reader for a FASTA file containing query data, parsing it into
/// [`QueryInfo`].
pub struct FastaQueryIter {
    reader: FastaReader<File>,
}

impl FastaQueryIter {
    /// Constructs a [`FastaQueryIter`] from an existing [`BufReader`].
    ///
    /// ## Errors
    ///
    /// See [`FastaReader::from_bufreader`].
    fn from_bufreader(buf: BufReader<File>) -> std::io::Result<Self> {
        let reader = FastaReader::from_bufreader(buf)?;
        Ok(Self { reader })
    }
}

impl Iterator for FastaQueryIter {
    type Item = std::io::Result<QueryInfo>;

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        let info = unwrap_or_return_some_err!(self.reader.next()?);
        let FastaSeq { name, sequence } = info;

        if name.contains('|') {
            let mut parts = name.split('|').map(|part| part.trim_ascii().to_string());

            if let (Some(id), Some(ctype)) = (parts.next(), parts.next())
                && !ctype.is_empty()
            {
                Some(Ok(QueryInfo {
                    id,
                    sequence,
                    ctype: Some(ctype),
                }))
            } else {
                Some(Err(std::io::Error::other(format!(
                    "Invalid FASTA header found. Expected ID or ID|ctype, found {name}"
                ))))
            }
        } else {
            Some(Ok(QueryInfo {
                id: name,
                sequence,
                ctype: None,
            }))
        }
    }
}

/// A reader for a TSV file containing query data, parsing it into
/// [`QueryInfo`].
///
/// Supports 3-column annotated (`ID<TAB>ctype<TAB>sequence`) and 2-column
/// unannotated (`ID\<TAB>sequence`) input.
pub struct TsvQueryIter {
    reader: Lines<BufReader<File>>,
}

impl TsvQueryIter {
    /// Constructs a [`TsvQueryIter`] from an already-initialized [`BufReader`].
    fn from_bufreader(buf: BufReader<File>) -> Self {
        Self { reader: buf.lines() }
    }

    /// Parses a single TSV line into a [`QueryInfo`].
    ///
    /// ## Validity
    ///
    /// The line should already have the `\n` or `\r\n` removed from the end,
    /// and should be non-empty and not solely contain whitespace.
    ///
    /// ## Errors
    ///
    /// The line must contain all required fields.
    fn parse_line(line: &str) -> std::io::Result<QueryInfo> {
        let mut columns = line.split('\t').map(|part| part.trim_ascii().to_string());

        // Validity: this will always exist since split is never empty
        let id = columns.next().unwrap_or_default();
        let Some(second) = columns.next() else {
            return Err(std::io::Error::other(
                "Invalid TSV format: expected 2 or 3 tab-separated columns, but found 1",
            ));
        };
        let third = columns.next();

        match third {
            // Three columns: ID  ctype  sequence  (annotated)
            Some(sequence) => {
                if second.is_empty() {
                    return Err(std::io::Error::other("Invalid TSV format: the second field was empty"));
                }

                let sequence = sequence.as_bytes().to_vec();

                Ok(QueryInfo {
                    id,
                    sequence,
                    ctype: Some(second),
                })
            }
            // Two columns: ID  sequence  (unannotated)
            None => {
                let sequence = second.as_bytes().to_vec();

                Ok(QueryInfo {
                    id,
                    sequence,
                    ctype: None,
                })
            }
        }
    }
}

impl Iterator for TsvQueryIter {
    type Item = std::io::Result<QueryInfo>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            let line = unwrap_or_return_some_err!(self.reader.next()?);

            if !line.trim().is_empty() {
                // Validity: Lines removes trailing line breaks, and we
                // ensure it is non-empty and not solely whitespace
                return Some(Self::parse_line(&line));
            }
        }
    }
}

define_whichever! {
    /// Unified query iterator over FASTA or TSV input, returning queries as
    /// [`QueryInfo`].
    pub enum QueryReader {
        /// Backed by a FASTA reader.
        Fasta(FastaQueryIter),
        /// Backed by a TSV reader.
        Tsv(TsvQueryIter),
    }

    impl Iterator for QueryReader {
        type Item = std::io::Result<QueryInfo>;
    }
}

impl QueryReader {
    /// Open `path`, peek at the first byte to detect format, and return the
    /// appropriate reader.
    ///
    /// ## Errors
    ///
    /// IO errors are propagated without context. An error is also returned if
    /// the file is empty, since the file is streamed in DAIS-ribosome and this
    /// error would not otherwise be detected.
    pub fn from_path(path: &Path) -> std::io::Result<Self> {
        let file = File::open(path).with_path_context("Failed to open file", path)?;

        let mut buffer = BufReader::new(file);

        match *buffer.peek(1)? {
            [] => Err(std::io::Error::other("The file is empty")),
            [b'>', ..] => Ok(QueryReader::Fasta(FastaQueryIter::from_bufreader(buffer)?)),
            _ => Ok(QueryReader::Tsv(TsvQueryIter::from_bufreader(buffer))),
        }
    }
}
