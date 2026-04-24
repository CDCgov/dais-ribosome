//! A [`Writers`] struct along with functions for writing all the TSV records
//! from DAIS-ribosome output.

use zoe::data::err::ResultWithErrorContext;

use crate::{
    data::ranges::StateRange,
    outputs::RibosomeOutput,
    tsv::{DelRowView, GenDelRowView, GenInsRowView, GenSeqRowView, InsRowView, SeqRowView},
};
use std::{
    fs::File,
    io::{BufWriter, Write},
    path::Path,
};

/// A set of three writers for sequence data, insertion data, and deletion data.
///
/// This same struct can be used for the product writers as well as the genome
/// writers.
#[derive(Clone, Eq, PartialEq, Hash, Debug)]
pub struct Writers<W> {
    /// The writer for the sequence output.
    pub seq: W,
    /// The writer for the insertion output.
    pub ins: W,
    /// The writer for the deletion output.
    pub del: W,
}

impl Writers<BufWriter<File>> {
    /// Constructs [`Writers`] by creating files at the three provided paths.
    ///
    /// ## Errors
    ///
    /// IO errors are propagated, with context including which argument had the
    /// error and the path.
    pub fn from_paths<P>(seq: P, ins: P, del: P) -> std::io::Result<Self>
    where
        P: AsRef<Path>, {
        let seq = seq.as_ref();
        let ins = ins.as_ref();
        let del = del.as_ref();

        let seq = File::create(seq).with_path_context("Could not open writer for sequence file", seq)?;
        let ins = File::create(ins).with_path_context("Could not open writer for insertion file", ins)?;
        let del = File::create(del).with_path_context("Could not open writer for deletion file", del)?;

        Ok(Self {
            seq: BufWriter::new(seq),
            ins: BufWriter::new(ins),
            del: BufWriter::new(del),
        })
    }
}

impl<W: Write> Writers<W> {
    /// Writes the product outputs for a single query to the appropriate writers.
    pub fn write_product_output(&mut self, output: &RibosomeOutput<'_>) -> std::io::Result<()> {
        for state in &output.states {
            for product in &state.products {
                let computed_product = &product.materialize(&output.query.nucleotides);

                let seq_row = SeqRowView::new(
                    computed_product,
                    &output.query.id,
                    &output.query.ctype,
                    state.reference_id,
                    output.formatting,
                );

                writeln!(self.seq, "{seq_row}")?;

                for insertion in &computed_product.insertions {
                    let ins_row = InsRowView::new(
                        insertion,
                        computed_product,
                        &output.query.id,
                        &output.query.ctype,
                        state.reference_id,
                    );

                    writeln!(self.ins, "{ins_row}")?;
                }

                for deletion in &computed_product.deletions {
                    let del_row = DelRowView::new(
                        deletion,
                        computed_product,
                        &output.query.id,
                        &output.query.ctype,
                        state.reference_id,
                    );

                    writeln!(self.del, "{del_row}")?;
                }
            }
        }

        Ok(())
    }

    /// Writes the genome outputs for a single query to the appropriate writers.
    pub fn write_genome_output(&mut self, output: &RibosomeOutput<'_>) -> std::io::Result<()> {
        for state in &output.states {
            let genome = state.materialize_genome(&output.query.nucleotides);

            let seq_row = GenSeqRowView::new(
                &genome,
                &output.query.id,
                &output.query.ctype,
                state.reference_id,
                output.formatting,
            );

            writeln!(self.seq, "{seq_row}")?;

            for insertion in &genome.insertions {
                let ins_row = GenInsRowView::new(insertion, &output.query.id, &output.query.ctype, state.reference_id);

                writeln!(self.ins, "{ins_row}")?;
            }

            for del_range in &state.genome_aln_states {
                if let StateRange::D(del) = del_range {
                    let del_row = GenDelRowView::new(del, &output.query.id, &output.query.ctype, state.reference_id);

                    writeln!(self.del, "{del_row}")?;
                }
            }
        }

        Ok(())
    }

    /// Transforms all the writers using a closure.
    pub fn map<U, F>(self, f: F) -> Writers<U>
    where
        F: Fn(W) -> U, {
        Writers {
            seq: f(self.seq),
            ins: f(self.ins),
            del: f(self.del),
        }
    }

    /// Flushes all the writers.
    ///
    /// This is also automatically called by [`write_product_output`] and
    /// [`write_genome_output`].
    ///
    /// [`write_product_output`]: Writers::write_product_output
    /// [`write_genome_output`]: Writers::write_genome_output
    pub fn flush(&mut self) -> std::io::Result<()> {
        self.seq.flush()?;
        self.ins.flush()?;
        self.del.flush()
    }
}
