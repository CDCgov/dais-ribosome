//! Configuration structs and specifications for DAIS-ribosome.

use crate::data::{exons::Exons, weights::CodonPositionWeights};

pub mod annotation_module;
pub(crate) mod cds_spec;
mod references;
pub mod toml;

/// The specifications for a single protein product (e.g., `HA`, `HA-signal`).
#[derive(Clone, Debug)]
pub(crate) struct ProductSpec {
    /// The protein/peptide product name
    pub(crate) name:          String,
    /// The exon coordinates, as well as translation rules specific to the exons
    pub(crate) exons:         Exons,
    /// The codon position weight matrix for the protein product
    pub(crate) codon_weights: Option<CodonPositionWeights>,
}
