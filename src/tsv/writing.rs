//! A [`Writers`] struct along with functions for writing all the TSV records
//! from DAIS-ribosome output.

use crate::{
    data::{RibosomeOutput, ranges::StateRange},
    error::RibosomeError,
    tsv::{DelRowView, GenDelRowView, GenInsRowView, GenSeqRowView, InsRowView, SeqRowView},
};
use std::io::Write;

/// A set of three writers for sequence data, insertion data, and deletion data.
///
/// This same struct can be used for the product writers as well as the genome
/// writers.
pub struct Writers<W> {
    /// The writer for the sequence output.
    pub seq: W,
    /// The writer for the insertion output.
    pub ins: W,
    /// The writer for the deletion output.
    pub del: W,
}

/// Writes the product outputs for a single query to the appropriate writers.
pub fn write_product_output<W: Write>(output: &RibosomeOutput<'_>, writers: &mut Writers<W>) -> Result<(), RibosomeError> {
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

            writeln!(writers.seq, "{seq_row}")?;

            for insertion in &computed_product.insertions {
                let ins_row = InsRowView::new(
                    insertion,
                    computed_product,
                    &output.query.id,
                    &output.query.ctype,
                    state.reference_id,
                );

                writeln!(writers.ins, "{ins_row}")?;
            }

            for deletion in &computed_product.deletions {
                let del_row = DelRowView::new(
                    deletion,
                    computed_product,
                    &output.query.id,
                    &output.query.ctype,
                    state.reference_id,
                );

                writeln!(writers.del, "{del_row}")?;
            }
        }
    }

    Ok(())
}

/// Writes the genome outputs for a single query to the appropriate writers.
pub fn write_genome_output<W: Write>(
    output: &RibosomeOutput<'_>, gen_writers: &mut Writers<W>,
) -> Result<(), RibosomeError> {
    for state in &output.states {
        let genome = state.materialize_genome(&output.query.nucleotides);

        let seq_row = GenSeqRowView::new(
            &genome,
            &output.query.id,
            &output.query.ctype,
            state.reference_id,
            output.formatting,
        );

        writeln!(gen_writers.seq, "{seq_row}")?;

        for insertion in &genome.insertions {
            let ins_row = GenInsRowView::new(insertion, &output.query.id, &output.query.ctype, state.reference_id);

            writeln!(gen_writers.ins, "{ins_row}")?;
        }

        for del_range in &state.genome_aln_states {
            if let StateRange::D(del) = del_range {
                let del_row = GenDelRowView::new(del, &output.query.id, &output.query.ctype, state.reference_id);

                writeln!(gen_writers.del, "{del_row}")?;
            }
        }
    }

    Ok(())
}
