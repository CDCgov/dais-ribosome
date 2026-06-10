//! A [`Writers`] struct along with functions for writing all the TSV records
//! from DAIS-ribosome output.

use crate::{
    data::ranges::StateRange,
    outputs::{MaybeComputedProduct, RibosomeOutput},
    tsv::{
        DelRowView, DelWriter, GenDelRowView, GenDelWriter, GenInsRowView, GenInsWriter, GenSeqRowView, GenSeqWriter,
        InsRowView, InsWriter, SeqRowView, SeqWriter,
    },
};
use std::{
    fs::File,
    io::{BufWriter, Write},
    path::Path,
};
use zoe::data::err::ResultWithErrorContext;

// TODO: Can this eventually use IRMA-core's trait? Right now IRMA-core is not a
// dependency
/// A trait allowing a writer to be finished, which may include flushing,
/// writing any footers, etc.
pub trait Finish {
    /// Finalizes the writer, performing flushing, writing any footers, etc.
    fn finish(self) -> std::io::Result<()>;
}

impl<W: Write> Finish for BufWriter<W> {
    fn finish(mut self) -> std::io::Result<()> {
        self.flush()
    }
}

impl<S, I, D> Finish for Writers<S, I, D>
where
    S: Finish,
    I: Finish,
    D: Finish,
{
    fn finish(self) -> std::io::Result<()> {
        self.seq.finish()?;
        self.ins.finish()?;
        self.del.finish()
    }
}

/// A convenience trait for a writer that is compatible with any of the output
/// records.
pub trait AnyWriter: SeqWriter + InsWriter + DelWriter + GenSeqWriter + GenInsWriter + GenDelWriter {}
impl<T: SeqWriter + InsWriter + DelWriter + GenSeqWriter + GenInsWriter + GenDelWriter> AnyWriter for T {}

/// A set of three writers for sequence data, insertion data, and deletion data.
///
/// This same struct can be used for the product writers as well as the genome
/// writers.
#[derive(Clone, Eq, PartialEq, Hash, Debug)]
pub struct Writers<S, I = S, D = S> {
    /// The writer for the sequence output.
    pub seq: S,
    /// The writer for the insertion output.
    pub ins: I,
    /// The writer for the deletion output.
    pub del: D,
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

impl<W> Writers<W> {
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
}

impl<S, I, D> Writers<S, I, D>
where
    S: SeqWriter,
    I: InsWriter,
    D: DelWriter,
{
    /// Writes the product outputs for a single query to the appropriate
    /// writers.
    pub fn write_product_output(&mut self, output: &RibosomeOutput<'_>) -> std::io::Result<()> {
        for state in &output.states {
            for product in &state.products {
                let computed_product = &product.materialize(&output.query);

                let seq_row = SeqRowView::new(
                    computed_product,
                    output.query.id(),
                    output.query.ctype(),
                    state.reference_id,
                    output.formatting,
                );

                match seq_row {
                    SeqRowView::Data(seq_data) => self.seq.write_seq_data(&seq_data)?,
                    SeqRowView::Empty(_) => {
                        // TODO: When we add the ability to toggle on null
                        // records, a method call will be needed here
                    }
                    SeqRowView::Deleted(deleted_seq_row) => self.seq.write_deleted_seq_row(&deleted_seq_row)?,
                }

                let computed_product = match computed_product {
                    MaybeComputedProduct::Ok(computed_product) => computed_product,
                    MaybeComputedProduct::Empty(_) => continue,
                    MaybeComputedProduct::Deleted(deleted_product) => {
                        let del_row = DelRowView::from_deleted_product(
                            deleted_product,
                            output.query.id(),
                            output.query.ctype(),
                            state.reference_id,
                        );

                        self.del.write_del_row(&del_row)?;
                        continue;
                    }
                };

                for insertion in &computed_product.insertions {
                    let ins_row = InsRowView::new(
                        insertion,
                        computed_product,
                        output.query.id(),
                        output.query.ctype(),
                        state.reference_id,
                    );

                    self.ins.write_ins_row(&ins_row)?;
                }

                for deletion in &computed_product.deletions {
                    let del_row = DelRowView::new(
                        deletion,
                        computed_product,
                        output.query.id(),
                        output.query.ctype(),
                        state.reference_id,
                    );

                    self.del.write_del_row(&del_row)?;
                }
            }
        }

        Ok(())
    }
}

impl<S, I, D> Writers<S, I, D>
where
    S: GenSeqWriter,
    I: GenInsWriter,
    D: GenDelWriter,
{
    /// Writes the genome outputs for a single query to the appropriate writers.
    pub fn write_genome_output(&mut self, output: &RibosomeOutput<'_>) -> std::io::Result<()> {
        for state in &output.states {
            let genome = state.materialize_genome(&output.query);

            let seq_row = GenSeqRowView::new(
                &genome,
                output.query.id(),
                output.query.ctype(),
                state.reference_id,
                output.formatting,
            );

            self.seq.write_gen_seq_row(&seq_row)?;

            for insertion in &genome.insertions {
                let ins_row = GenInsRowView::new(insertion, output.query.id(), output.query.ctype(), state.reference_id);

                self.ins.write_gen_ins_row(&ins_row)?;
            }

            for del_range in &state.genome_aln_states {
                if let StateRange::D(del) = del_range {
                    let del_row = GenDelRowView::new(del, output.query.id(), output.query.ctype(), state.reference_id);

                    self.del.write_gen_del_row(&del_row)?;
                }
            }
        }

        Ok(())
    }
}
