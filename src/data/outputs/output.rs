use crate::{
    config::Formatting,
    data::{
        ComputedGenomeInsertion, DelRow, GenDelRow, GenInsRow, GenRow, InsRow, PrecomputedGenomeData, QueryRecord, SeqRow,
        products::{ComputedProduct, Product},
        ranges::{InsertionIdx, StateRange},
    },
    hashing::nt_id,
};
use std::{ops::Range, sync::OnceLock};
use zoe::prelude::*;

#[derive(Debug)]
pub struct RibosomeOutput<'a> {
    /// Original query record
    pub query:             QueryRecord,
    /// Both genome and protein product alignment states
    pub(crate) states:     Vec<GenomeAndProductStates<'a>>,
    /// Output formatting rules, parsed from the TOML configuration.
    pub(crate) formatting: &'a Formatting,
}

impl<'a> RibosomeOutput<'a> {
    pub fn materialize(self) -> ComputedRibosomeOutput<'a> {
        // Validity: query.nucleotides contains unaligned, uppercase IUPAC bases
        ComputedRibosomeOutput {
            states:     self
                .states
                .into_iter()
                .map(|state| state.materialize(&self.query.nucleotides))
                .collect(),
            query:      self.query,
            formatting: self.formatting,
        }
    }
}

#[derive(Debug)]
pub struct ComputedRibosomeOutput<'a> {
    /// Original query record
    pub query:             QueryRecord,
    /// Both genome and protein product alignment states
    pub(crate) states:     Vec<ComputedGenomeAndProductStates<'a>>,
    /// Output formatting rules
    pub(crate) formatting: &'a Formatting,
}

#[derive(Debug)]
pub(crate) struct GenomeAndProductStates<'a> {
    /// Reference ID
    pub(crate) reference_id: &'a str,

    /// The length of the reference sequence.
    pub(crate) ref_len: usize,

    /// Genome alignment to nucleotide reference sequence expressed as [`StateRange`]
    pub(crate) genome_aln_states: Vec<StateRange>,

    /// The range of the stop extension within the query, if present.
    pub(crate) stop_extension_query_range: Option<Range<usize>>,

    /// The number of bases in the reference sequence that were not aligned
    /// against in the beginning (i.e., not included in `genome_aln_states`).
    pub(crate) leading_ref_unaligned: usize,

    /// The number of bases in the reference sequence that were not aligned
    /// against at the end (i.e., not included in `genome_aln_states`).
    pub(crate) trailing_ref_unaligned: usize,

    /// Contains all relevant product data, including the protein name.
    pub(crate) products: Vec<Product<'a>>,

    /// Lazily computed genome data, cached via OnceLock.
    pub(crate) computed_genome: OnceLock<PrecomputedGenomeData>,
}

impl<'a> GenomeAndProductStates<'a> {
    /// Computes the output data for all products in TODO
    ///
    /// ## Validity
    ///
    /// The `query` should contain unaligned, uppercase IUPAC bases.
    fn materialize(self, query: &Nucleotides) -> ComputedGenomeAndProductStates<'a> {
        ComputedGenomeAndProductStates {
            reference_id:               self.reference_id,
            ref_len:                    self.ref_len,
            genome_aln_states:          self.genome_aln_states,
            stop_extension_query_range: self.stop_extension_query_range,
            leading_ref_unaligned:      self.leading_ref_unaligned,
            trailing_ref_unaligned:     self.trailing_ref_unaligned,
            computed_genome:            self.computed_genome,
            products:                   self.products.into_iter().map(|product| product.materialize(query)).collect(),
        }
    }
}

#[derive(Debug)]
pub(crate) struct ComputedGenomeAndProductStates<'a> {
    /// The reference ID.
    pub(crate) reference_id: &'a str,

    /// The length of the reference sequence.
    pub(crate) ref_len: usize,

    /// Genome alignment to nucleotide reference sequence expressed as
    /// [`StateRange`].
    pub(crate) genome_aln_states: Vec<StateRange>,

    /// The range of the stop extension within the query, if present.
    pub(crate) stop_extension_query_range: Option<Range<usize>>,

    /// The number of bases in the reference sequence that were not aligned
    /// against in the beginning (i.e., not included in `genome_aln_states`).
    pub(crate) leading_ref_unaligned: usize,

    /// The number of bases in the reference sequence that were not aligned
    /// against at the end (i.e., not included in `genome_aln_states`).
    pub(crate) trailing_ref_unaligned: usize,

    /// Contains all relevant product data, including the protein name.
    pub(crate) products: Vec<ComputedProduct<'a>>,

    /// Lazily computed genome data, cached via OnceLock.
    pub(crate) computed_genome: OnceLock<PrecomputedGenomeData>,
}

impl<'a> ComputedRibosomeOutput<'a> {
    /// Generates `.seq` output rows for this query.
    pub fn seq_rows(&self) -> impl Iterator<Item = SeqRow<'_>> + '_ {
        self.states.iter().flat_map(move |state| {
            state.products.iter().map(move |computed_product| SeqRow {
                id: &self.query.id,
                ctype: &self.query.ctype,
                ref_id: state.reference_id,
                computed_product,
                formatting: self.formatting,
            })
        })
    }

    /// Generates `.ins` output rows for this query.
    pub fn ins_rows(&self) -> impl Iterator<Item = InsRow<'_>> + '_ {
        self.states.iter().flat_map(move |state| {
            state.products.iter().flat_map(move |computed_product| {
                computed_product.insertions.iter().map(move |insertion| InsRow {
                    id: &self.query.id,
                    ctype: &self.query.ctype,
                    ref_id: state.reference_id,
                    protein: computed_product.product_name,
                    insertion,
                })
            })
        })
    }

    /// Generates displayable deletion rows
    pub fn del_rows<'b>(&'b self) -> impl Iterator<Item = DelRow<'b>> + 'b {
        self.states.iter().flat_map(move |state| {
            state.products.iter().flat_map(move |computed_product| {
                computed_product.deletions.iter().map(move |deletion| DelRow {
                    id: &self.query.id,
                    ctype: &self.query.ctype,
                    ref_id: state.reference_id,
                    protein: computed_product.product_name,
                    computed_product,
                    deletion,
                })
            })
        })
    }

    /// Eagerly materialize all products and genome data for use in
    /// parallel executors.
    pub fn materialize(&self) {
        for state in &self.states {
            state.materialize_genome(&self.query.nucleotides);
        }
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

impl<'a> ComputedGenomeAndProductStates<'a> {
    /// Lazily compute and cache genome data from genome alignment states.
    pub fn materialize_genome(&self, query: &Nucleotides) -> &PrecomputedGenomeData {
        self.computed_genome.get_or_init(|| {
            let mut genome_seq = Nucleotides::new();
            let mut genome_aln = Nucleotides::from(vec![b'.'; self.leading_ref_unaligned]);
            let mut insertions = Vec::new();
            let mut has_insertion = false;

            for state in &self.genome_aln_states {
                match state {
                    StateRange::M(m) => {
                        let slice = &query[m.query_range.clone()];
                        genome_seq.extend_from_slice(slice);
                        genome_aln.extend_from_slice(slice);
                    }
                    StateRange::I(ins) => {
                        let slice = &query[ins.query_range.clone()];
                        genome_seq.extend_from_slice(slice);
                        has_insertion = true;
                        insertions.push(ComputedGenomeInsertion::new(ins.ref_index, slice));
                    }
                    StateRange::D(del) => {
                        genome_aln.pad_end(b'-', del.ref_range.len());
                    }
                }
            }

            if let Some(ref ext_range) = self.stop_extension_query_range {
                let nt_insertion_idx = InsertionIdx::from_right_idx(self.ref_len);
                let slice = &query[ext_range.clone()];
                genome_seq.extend_from_slice(slice);
                has_insertion = true;
                insertions.push(ComputedGenomeInsertion::new(nt_insertion_idx, slice));
            }

            let genome_id = nt_id(&genome_seq);
            let genome_length = genome_seq.len();

            PrecomputedGenomeData {
                genome_id,
                genome_length,
                has_insertion,
                genome_seq,
                genome_aln,
                insertions,
                trailing_ref_unaligned: self.trailing_ref_unaligned,
            }
        })
    }
}
