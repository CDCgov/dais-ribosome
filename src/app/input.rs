//! Multiplexed query input: auto-detects FASTA vs. TSV by peeking at the
//! first byte, then yields [`QueryRecord`]s through a unified iterator.
//! The file is opened only once.

use dais_ribosome::{QueryRecord, error::RibosomeError};
use std::{
    fs::File,
    io::{BufRead, BufReader, Lines},
    path::Path,
};
use zoe::{
    data::{ByteMap, RetainSequence, err::ResultWithErrorContext, fasta::FastaSeq},
    define_whichever,
    prelude::*,
};

use crate::app::log::time_stamp;

/// A reader for a FASTA file containing query data, parsing it into
/// [`QueryRecord`].
///
/// This uses *Zoe*'s [`FastaReader`]. All non-IUPAC bases and gap characters
/// are filtered, and bases are converted to uppercase. `U` is preserved in
/// addition to `T`.
pub struct FastaQueryIter {
    reader: FastaReader<File>,
}

impl FastaQueryIter {
    fn from_bufreader(buf: BufReader<File>) -> Result<Self, RibosomeError> {
        let reader = FastaReader::from_bufreader(buf)?;
        Ok(Self { reader })
    }
}

impl Iterator for FastaQueryIter {
    type Item = Result<QueryRecord, QueryInputError>;

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        let record = match self.reader.next()? {
            Ok(record) => record,
            Err(e) => return Some(Err(e.into())),
        };

        Some(parse_fasta_seq(record))
    }
}

pub enum QueryInputError {
    Io(std::io::Error),
    NoNucleotides(String, ReaderType),
}

pub enum ReaderType {
    Fasta,
    Tsv,
}

impl From<std::io::Error> for QueryInputError {
    fn from(value: std::io::Error) -> Self {
        QueryInputError::Io(value)
    }
}

impl From<&str> for QueryInputError {
    fn from(value: &str) -> Self {
        QueryInputError::Io(std::io::Error::other(value))
    }
}

impl From<String> for QueryInputError {
    fn from(value: String) -> Self {
        QueryInputError::Io(std::io::Error::other(value))
    }
}

/// Parses a [`FastaSeq`] into a [`QueryRecord`] for use in DAIS-ribosome.
///
///
/// This filters all non-IUPAC bases and gap characters, and converts to
/// uppercase. `U` is preserved in addition to `T`.
///
/// ## Errors
///
/// - The sequence must be non-empty after filtering.
/// - The header must be successfully parsed.
fn parse_fasta_seq(record: FastaSeq) -> Result<QueryRecord, QueryInputError> {
    let FastaSeq { name, sequence } = record;
    // BREAKING: we previously only removed: '*: .~-'
    let nucleotides = sanitize_seq(sequence);

    if nucleotides.is_empty() {
        return Err(QueryInputError::NoNucleotides(name, ReaderType::Fasta));
    }

    if name.contains('|') {
        let mut parts = name.split('|').map(|part| part.trim_ascii().to_string());

        if let (Some(id), Some(ctype)) = (parts.next(), parts.next())
            && !ctype.is_empty()
        {
            Ok(QueryRecord { id, nucleotides, ctype })
        } else {
            Err(format!("Invalid FASTA header found. Expected ID or ID|ctype, found {name}").into())
        }
    } else {
        // TODO : handle unclassified queries
        Ok(QueryRecord {
            id: name.trim_ascii().to_string(),
            nucleotides,
            ctype: String::new(),
        })
    }
}

/// A reader for a TSV file containing query data, parsing it into
/// [`QueryRecord`].
///
/// Supports 3-column annotated (`ID<TAB>ctype<TAB>sequence`) and 2-column
/// unannotated (`ID\<TAB>sequence`, currently stubbed) input.
///
/// All non-IUPAC bases and gap characters are filtered, and bases are converted
/// to uppercase. `U` is preserved in addition to `T`.
pub struct TsvQueryIter {
    reader: Lines<BufReader<File>>,
}

impl TsvQueryIter {
    /// Constructs a [`TsvQueryIter`] from an already-initialized [`BufReader`].
    fn from_bufreader(buf: BufReader<File>) -> Self {
        Self { reader: buf.lines() }
    }

    /// Parses a single TSV line into a [`QueryRecord`].
    ///
    /// This filters all non-IUPAC bases and gap characters, and converts to
    /// uppercase. `U` is preserved in addition to `T`.
    ///
    /// ## Validity
    ///
    /// The line should already have the `\n` or `\r\n` removed from the end,
    /// and should be non-empty and not solely contain whitespace.
    ///
    /// ## Errors
    ///
    /// If the line is missing fields, or if the sequence contained no unaligned
    /// DNA data, a [`RibosomeError::Io`] error is returned.
    fn parse_line(line: &str) -> Result<QueryRecord, QueryInputError> {
        let mut columns = line.split('\t');

        // Validity: this will always exist since split is never empty
        let id = columns.next().unwrap_or_default();
        let Some(second) = columns.next() else {
            return Err("Invalid TSV format: expected 2 or 3 tab-separated columns, but found 1".into());
        };
        let third = columns.next();

        match third {
            // Three columns: ID  ctype  sequence  (annotated)
            Some(seq_field) => {
                let ctype = second.trim_ascii().to_string();
                if ctype.is_empty() {
                    return Err("Invalid TSV format: the second field was empty".into());
                }

                let nucleotides = sanitize_seq(seq_field.as_bytes().to_vec());

                if nucleotides.is_empty() {
                    return Err(QueryInputError::NoNucleotides(id.to_string(), ReaderType::Tsv));
                }

                Ok(QueryRecord {
                    id: id.trim_ascii().to_string(),
                    nucleotides,
                    ctype,
                })
            }
            // Two columns: ID  sequence  (unannotated — stub)
            None => {
                let nucleotides = sanitize_seq(second.as_bytes().to_vec());

                if nucleotides.is_empty() {
                    return Err(QueryInputError::NoNucleotides(id.to_string(), ReaderType::Tsv));
                }

                // TODO: handle unclassified TSV queries (feature-gated)
                Ok(QueryRecord {
                    id: id.trim_ascii().to_string(),
                    nucleotides,
                    ctype: String::new(),
                })
            }
        }
    }
}

impl Iterator for TsvQueryIter {
    type Item = Result<QueryRecord, QueryInputError>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            let line = match self.reader.next()? {
                Ok(line) => line,
                Err(e) => return Some(Err(e.into())),
            };

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
    /// [`QueryRecord`].
    ///
    /// Both iterators filter all non-IUPAC bases and gap characters, and bases
    /// are converted to uppercase. `U` is preserved in addition to `T`.
    pub enum QueryInput {
        /// Backed by a FASTA reader.
        Fasta(FastaQueryIter),
        /// Backed by a TSV reader.
        Tsv(TsvQueryIter),
    }

    impl Iterator for QueryInput {
        type Item = Result<QueryRecord, QueryInputError>;
    }
}

impl QueryInput {
    /// Open `path`, peek at the first byte to detect format, and return
    /// the appropriate reader.
    pub fn open(path: &Path) -> Result<Self, RibosomeError> {
        let file = File::open(path).with_path_context("Failed to open file", path)?;

        let mut buffer = BufReader::new(file);

        match *buffer.peek(1)? {
            [] => Err(RibosomeError::EmptyFile(path.to_path_buf())),
            [b'>', ..] => Ok(QueryInput::Fasta(FastaQueryIter::from_bufreader(buffer)?)),
            _ => Ok(QueryInput::Tsv(TsvQueryIter::from_bufreader(buffer))),
        }
    }
}

/// An extension trait for handling [`QueryInputError::NoNucleotides`] in an
/// iterator.
pub trait HandleNoNucleotidesExt {
    /// Filters any records with empty sequences post-filtering, issuing a
    /// warning to `stderr` if `warn_no_nucleotides` is set.
    fn handle_no_nucleotides(self, warn_no_nucleotides: bool) -> impl Iterator<Item = std::io::Result<QueryRecord>>;
}

impl<I> HandleNoNucleotidesExt for I
where
    I: Iterator<Item = Result<QueryRecord, QueryInputError>>,
{
    fn handle_no_nucleotides(self, warn_no_nucleotides: bool) -> impl Iterator<Item = std::io::Result<QueryRecord>> {
        self.filter_map(move |res| match res {
            Ok(record) => Some(Ok(record)),
            Err(QueryInputError::Io(e)) => Some(Err(e)),
            Err(QueryInputError::NoNucleotides(header, reader_type)) => {
                if warn_no_nucleotides {
                    let field = match reader_type {
                        ReaderType::Fasta => "header",
                        ReaderType::Tsv => "ID",
                    };

                    time_stamp(
                        &format!("A sequence contained no unaligned DNA data. See {field}: {header}"),
                        true,
                    );
                }
                None
            }
        })
    }
}

/// Sanitizes an incoming sequence so that it meets the validity requirements of
/// [`QueryRecord`].
///
/// This converts to uppercase, preserves IUPAC characters, preserves `U` in
/// addition to `T`, and preserves `X`. All other bytes are removed.
#[must_use]
#[cfg(feature = "regression-testing")]
fn sanitize_seq(mut seq: Vec<u8>) -> Nucleotides {
    const SANITIZE: ByteMap = ByteMap::all(0)
        .preserve_range(b'A'..=b'Z')
        .preserve_range(b'a'..=b'z')
        .map(b"acgturyswkmbdhvn", b"ACGTURYSWKMBDHVN")
        .map(b"ux", b"UX");

    seq.retain_by_recoding(&SANITIZE);
    Nucleotides::from(seq)
}

/// Sanitizes an incoming sequence so that it meets the validity requirements of
/// [`QueryRecord`].
///
/// This converts to uppercase, preserves IUPAC characters, preserves `U` in
/// addition to `T`, and maps `X` to `N`. All other bytes are removed.
#[must_use]
#[cfg(not(feature = "regression-testing"))]
fn sanitize_seq(mut seq: Vec<u8>) -> Nucleotides {
    const SANITIZE: ByteMap = ByteMap::all(0)
        .preserve_range(b'A'..=b'Z')
        .preserve_range(b'a'..=b'z')
        .map(b"acgturyswkmbdhvn", b"ACGTURYSWKMBDHVN")
        .map(b"uxX", b"UNN");

    seq.retain_by_recoding(&SANITIZE);
    Nucleotides::from(seq)
}
