//! Multiplexed query input: auto-detects FASTA vs. TSV by peeking at the
//! first byte, then yields [`QueryRecord`]s through a unified iterator.
//! The file is opened only once.

use dais_ribosome::data::{QueryRecord, RibosomeError};
use std::{
    fs::File,
    io::{BufRead, BufReader},
    path::Path,
};
use zoe::{data::nucleotides::ToDNA, prelude::*};

/// FASTA → [`QueryRecord`] adapter over Zoe's [`FastaReader`].
pub struct FastaQueryIter {
    inner: FastaReader<File>,
}

impl FastaQueryIter {
    fn from_bufreader(buf: BufReader<File>) -> Result<Self, RibosomeError> {
        let inner = FastaReader::from_bufreader(buf)?;
        Ok(Self { inner })
    }
}

impl Iterator for FastaQueryIter {
    type Item = Result<QueryRecord, RibosomeError>;

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        let fasta_seq = self.inner.next()?;
        Some(fasta_seq.map_err(RibosomeError::from).and_then(QueryRecord::try_from))
    }
}

/// TSV → [`QueryRecord`] line reader.
///
/// Supports 3-column annotated (`ID\tctype\tsequence`) and 2-column
/// unannotated (`ID\tsequence`, currently stubbed) input.
pub struct TsvQueryIter {
    reader: BufReader<File>,
    line:   String,
}

impl TsvQueryIter {
    fn from_bufreader(buf: BufReader<File>) -> Self {
        Self {
            reader: buf,
            line:   String::new(),
        }
    }

    /// Parse a single TSV line into a [`QueryRecord`].
    fn parse_line(line: &str) -> Result<QueryRecord, RibosomeError> {
        let line = line.trim_end_matches(['\r', '\n']);

        if line.is_empty() || line.bytes().all(|b| b.is_ascii_whitespace()) {
            return Err(RibosomeError::InvalidTsvFormat);
        }

        let mut columns = line.split('\t');

        let id = columns.next().unwrap_or_default(); // always exists after split
        let second = columns.next().ok_or(RibosomeError::InvalidTsvFormat)?;
        let third = columns.next();

        match third {
            // Three columns: ID  ctype  sequence  (annotated)
            Some(seq_field) => {
                let ctype = second.trim_ascii().to_string();
                if ctype.is_empty() {
                    return Err(RibosomeError::NoCtype(id.to_string()));
                }

                let nucleotides = seq_field.as_bytes().to_vec().filter_to_dna_unaligned();
                if nucleotides.is_empty() {
                    return Err(RibosomeError::InvalidSequence(id.to_string()));
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
                    return Err(RibosomeError::InvalidSequence(id.to_string()));
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
    type Item = Result<QueryRecord, RibosomeError>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            self.line.clear();
            match self.reader.read_line(&mut self.line) {
                Ok(0) => return None,
                Ok(_) => {
                    // Skip blank lines
                    if self.line.trim().is_empty() {
                        continue;
                    }
                    return Some(Self::parse_line(&self.line));
                }
                Err(e) => return Some(Err(RibosomeError::IO(e))),
            }
        }
    }
}

zoe::define_whichever! {
    /// Unified query iterator over FASTA or TSV input.
    pub enum QueryInput {
        /// Backed by a FASTA reader.
        Fasta(FastaQueryIter),
        /// Backed by a TSV reader.
        Tsv(TsvQueryIter),
    }

    impl Iterator for QueryInput {
        type Item = Result<QueryRecord, RibosomeError>;
    }
}

impl QueryInput {
    /// Open `path`, peek at the first byte to detect format, and return
    /// the appropriate reader.
    pub fn open(path: &Path) -> Result<Self, RibosomeError> {
        // TODO: Check non-empty
        let file = File::open(path)?;
        let mut buffer = BufReader::new(file);

        let first = buffer.peek(1)?;
        if first.is_empty() {
            return Err(RibosomeError::BlankFirstLine(path.to_path_buf()));
        }

        match first[0] {
            b'>' => Ok(QueryInput::Fasta(FastaQueryIter::from_bufreader(buffer)?)),
            _ => Ok(QueryInput::Tsv(TsvQueryIter::from_bufreader(buffer))),
        }
    }
}
