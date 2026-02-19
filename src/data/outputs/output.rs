use crate::{
    annotation::hashing::nt_id,
    config::Formatting,
    data::{
        ComputedGenomeInsertion, DelRow, GenDelRow, GenInsRow, GenRow, InsRow, PrecomputedGenomeData, QueryRecord, SeqRow,
        products::Product, ranges::StateRange,
    },
};
use std::sync::OnceLock;
use zoe::prelude::*;

#[derive(Debug)]
pub struct RibosomeOutput<'a> {
    /// Original query record
    pub query:             QueryRecord,
    /// Both genome and protein product alignment states
    pub(crate) states:     Vec<GenomeAndProductStates<'a>>,
    /// Output formatting rules
    pub(crate) formatting: &'a Formatting,
}

#[derive(Debug)]
pub(crate) struct GenomeAndProductStates<'a> {
    /// Reference ID
    pub(crate) reference_id:      &'a str,
    /// Reference sequence length for genome padding
    pub(crate) ref_len:           usize,
    /// Genome alignment to nucleotide reference sequence expressed as [`StateRange`]
    pub(crate) genome_aln_states: Vec<StateRange>,
    /// Contains all relevant product data, including the protein name.
    pub(crate) products:          Vec<Product<'a>>,
    /// Lazily computed genome data, cached via OnceLock.
    pub(crate) computed_genome:   OnceLock<PrecomputedGenomeData>,
}

impl<'a> RibosomeOutput<'a> {
    /// Eagerly materialize all products and genome data for use in
    /// parallel executors.
    pub fn materialize(&self) {
        for state in &self.states {
            state.materialize_genome(&self.query.nucleotides);
            for product in &state.products {
                product.materialize(&self.query.nucleotides);
            }
        }
    }

    /// Generate `.seq` output rows for this query.
    pub fn seq_rows(&self) -> impl Iterator<Item = SeqRow<'_>> + '_ {
        self.states.iter().flat_map(move |state| {
            state.products.iter().map(move |product| SeqRow {
                id: &self.query.id,
                ctype: &self.query.ctype,
                ref_id: state.reference_id,
                product,
                computed_product: product.materialize(&self.query.nucleotides),
                formatting: self.formatting,
            })
        })
    }

    /// Generate `.ins` output rows for this query.
    pub fn ins_rows(&self) -> impl Iterator<Item = InsRow<'_>> + '_ {
        self.states.iter().flat_map(move |state| {
            state.products.iter().flat_map(move |product| {
                let computed = product.materialize(&self.query.nucleotides);
                computed
                    .insertions
                    .iter()
                    .filter(|ins| !ins.filtered)
                    .map(move |insertion| InsRow {
                        id: &self.query.id,
                        ctype: &self.query.ctype,
                        ref_id: state.reference_id,
                        protein: &product.product_spec.name,
                        insertion,
                    })
            })
        })
    }

    /// Generate displayable deletion rows
    pub fn del_rows<'b>(&'b self) -> impl Iterator<Item = DelRow<'b>> + 'b {
        self.states.iter().flat_map(move |state| {
            state.products.iter().flat_map(move |product| {
                let computed = product.materialize(&self.query.nucleotides);
                computed.deletions.iter().map(move |deletion| DelRow {
                    id: &self.query.id,
                    ctype: &self.query.ctype,
                    ref_id: state.reference_id,
                    protein: &product.product_spec.name,
                    computed_product: computed,
                    deletion,
                })
            })
        })
    }

    /// Generate genome alignment rows for display
    pub fn gen_rows(&self) -> impl Iterator<Item = GenRow<'_>> + '_ {
        self.states.iter().map(move |state| {
            let genome = state.materialize_genome(&self.query.nucleotides);
            GenRow {
                id: &self.query.id,
                ctype: &self.query.ctype,
                ref_id: state.reference_id,
                genome,
                ref_len: state.ref_len,
                formatting: self.formatting,
            }
        })
    }

    /// Generate genome insertion rows for display
    pub fn gen_ins_rows<'b>(&'b self) -> impl Iterator<Item = GenInsRow<'b>> + 'b {
        self.states.iter().flat_map(move |state| {
            let genome = state.materialize_genome(&self.query.nucleotides);
            genome.insertions.iter().map(move |insertion| GenInsRow {
                id: &self.query.id,
                ctype: &self.query.ctype,
                ref_id: state.reference_id,
                insertion,
            })
        })
    }

    /// Generate genome deletion rows for display
    pub fn gen_del_rows<'b>(&'b self) -> impl Iterator<Item = GenDelRow<'b>> + 'b {
        self.states.iter().flat_map(move |state| {
            state.genome_aln_states.iter().filter_map(move |range| {
                if let StateRange::D(del) = range {
                    Some(GenDelRow {
                        id:       &self.query.id,
                        ctype:    &self.query.ctype,
                        ref_id:   state.reference_id,
                        deletion: del,
                    })
                } else {
                    None
                }
            })
        })
    }
}

impl<'a> GenomeAndProductStates<'a> {
    /// Lazily compute and cache genome data from genome alignment states.
    pub fn materialize_genome(&self, query: &Nucleotides) -> &PrecomputedGenomeData {
        self.computed_genome.get_or_init(|| {
            let mut genome_seq_bytes = Vec::new();
            let mut genome_aln_bytes = Vec::new();
            let mut insertions = Vec::new();
            let mut has_insertion = false;

            let query_bytes = query.as_bytes();

            for state in &self.genome_aln_states {
                match state {
                    StateRange::M(m) => {
                        let slice = &query_bytes[m.query_range.clone()];
                        genome_seq_bytes.extend_from_slice(slice);
                        genome_aln_bytes.extend_from_slice(slice);
                    }
                    StateRange::I(ins) => {
                        let slice = &query_bytes[ins.query_range.clone()];
                        genome_seq_bytes.extend_from_slice(slice);
                        has_insertion = true;

                        let inserted_nucleotides = Nucleotides::from_vec_unchecked(slice.to_vec());
                        insertions.push(ComputedGenomeInsertion {
                            upstream_nt: ins.upstream_ref_index + 1,
                            inserted_nucleotides,
                        });
                    }
                    StateRange::D(del) => {
                        genome_aln_bytes.extend(std::iter::repeat_n(b'-', del.ref_range.len()));
                    }
                }
            }

            let genome_seq = Nucleotides::from_vec_unchecked(genome_seq_bytes);
            let genome_id = nt_id(&genome_seq).unwrap_or_default();
            let genome_length = genome_seq.len();
            let genome_aln = Nucleotides::from_vec_unchecked(genome_aln_bytes);

            PrecomputedGenomeData {
                genome_id,
                genome_length,
                has_insertion,
                genome_seq,
                genome_aln,
                insertions,
            }
        })
    }
}
