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
    data::{err::ResultWithErrorContext, fasta::FastaSeq, nucleotides::ToDNA},
    define_whichever,
    prelude::*,
    unwrap_or_return_some_err,
};

/// A reader for a FASTA file containing query data, parsing it into
/// [`QueryRecord`].
///
/// This uses *Zoe*'s [`FastaReader`].
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
    type Item = std::io::Result<QueryRecord>;

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        Some(self.reader.next()?.and_then(parse_fasta_seq))
    }
}

fn parse_fasta_seq(record: FastaSeq) -> std::io::Result<QueryRecord> {
    let FastaSeq { name, sequence } = record;
    // BREAKING: we previously only removed: '*: .~-'
    let nucleotides = sequence.filter_to_dna_unaligned();

    if nucleotides.is_empty() {
        return Err(std::io::Error::other(format!(
            "A sequence contained no unaligned DNA data. See header: {name}"
        )));
    }

    if name.contains('|') {
        let mut parts = name.split('|').map(|part| part.trim_ascii().to_string());

        if let (Some(id), Some(ctype)) = (parts.next(), parts.next())
            && !ctype.is_empty()
        {
            Ok(QueryRecord { id, nucleotides, ctype })
        } else {
            Err(std::io::Error::other(format!(
                "Invalid FASTA header found. Expected ID or ID|ctype, found {name}"
            )))
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
    /// ## Validity
    ///
    /// The line should already have the `\n` or `\r\n` removed from the end,
    /// and should be non-empty and not solely contain whitespace.
    ///
    /// ## Errors
    ///
    /// If the line is missing fields, or if the sequence contained no unaligned
    /// DNA data, a [`RibosomeError::Io`] error is returned.
    fn parse_line(line: &str) -> std::io::Result<QueryRecord> {
        let mut columns = line.split('\t');

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
            Some(seq_field) => {
                let ctype = second.trim_ascii().to_string();
                if ctype.is_empty() {
                    return Err(std::io::Error::other("Invalid TSV format: the second field was empty"));
                }

                let nucleotides = seq_field.as_bytes().to_vec().filter_to_dna_unaligned();

                if nucleotides.is_empty() {
                    return Err(std::io::Error::other(format!(
                        "A sequence contained no unaligned DNA data. See ID: {id}"
                    )));
                }

                Ok(QueryRecord {
                    id: id.trim_ascii().to_string(),
                    nucleotides,
                    ctype,
                })
            }
            // Two columns: ID  sequence  (unannotated — stub)
            None => {
                let nucleotides = second.as_bytes().to_vec().filter_to_dna_unaligned();

                if nucleotides.is_empty() {
                    return Err(std::io::Error::other(format!(
                        "A sequence contained no unaligned DNA data. See ID: {id}"
                    )));
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
    type Item = std::io::Result<QueryRecord>;

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
    /// Unified query iterator over FASTA or TSV input.
    pub enum QueryInput {
        /// Backed by a FASTA reader.
        Fasta(FastaQueryIter),
        /// Backed by a TSV reader.
        Tsv(TsvQueryIter),
    }

    impl Iterator for QueryInput {
        type Item = std::io::Result<QueryRecord>;
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
