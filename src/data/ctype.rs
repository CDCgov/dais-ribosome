//! Compound type hierarchy: ctype → reference_id → protein_product.
use crate::{
    config::{AlignmentParams, product_spec::ProductSpec},
    data::{
        exons::Exons,
        keys::{RefKey, SpecKey},
        weights::CodonWeightMatrix,
    },
};
use std::collections::HashMap;
use zoe::{
    alignment::{Alignment, SharedProfiles},
    data::{err::ResultWithErrorContext, nucleotides::Nucleotides},
    iter_utils::ProcessResultsExt,
    prelude::*,
};

/// Top-level index: ctype string → compound type data.
pub type CompoundTypeMap<'a> = HashMap<String, Vec<ReferenceGroup<'a>>>;

/// Pre-computed alignment profile for a reference sequence.
pub type AlignmentProfiles<'a> = SharedProfiles<'a, 32, 16, 8, 5>;

/// Information about references sharing the same `reference_id` within a
/// compound type.
///
/// All the references must be the same length.
#[derive(Debug)]
pub(crate) struct ReferenceGroup<'a> {
    /// The shared reference ID of the reference sequences.
    pub(crate) reference_id: String,
    /// The shared length of the reference sequences.
    pub(crate) length:       usize,
    /// The alignment profiles corresponding to the sequences.
    pub(crate) profiles:     Vec<AlignmentProfiles<'a>>,
    pub(crate) proteins:     Vec<ProductSpec<'a>>,
}

impl<'a> ReferenceGroup<'a> {
    pub fn new(
        ref_key: &RefKey, seqs: &'a [Nucleotides], params: &'a AlignmentParams,
        cds_spec: &'a HashMap<RefKey, Vec<(String, Exons)>>, codon_weights: &'a CodonWeightMatrix,
    ) -> std::io::Result<Self> {
        let length = seqs.first().map_or(0, Nucleotides::len);
        for seq in seqs {
            if seq.len() != length {
                return Err(std::io::Error::other(format!(
                    "Inconsistent reference lengths for '{reference_id}|{compound_type}'",
                    reference_id = ref_key.reference_id,
                    compound_type = ref_key.compound_type
                )));
            }
        }

        let profiles = seqs
            .iter()
            .map(|seq| {
                AlignmentProfiles::new(seq, &params.matrix, params.gap_open, params.gap_extend)
                    .with_context("Failed to build alignment profiles")
            })
            .collect::<Result<Vec<_>, _>>()?;

        let proteins = cds_spec
            .get(ref_key)
            .into_iter()
            .flatten()
            .map(|(protein_name, exons)| {
                let spec_key = SpecKey::new(&ref_key.reference_id, protein_name);
                ProductSpec {
                    name: protein_name,
                    exons,
                    codon_weights: codon_weights.get(&spec_key),
                }
            })
            .collect();

        Ok(Self {
            reference_id: ref_key.reference_id.to_string(),
            length,
            profiles,
            proteins,
        })
    }

    pub fn extend(&mut self, seqs: &'a [Nucleotides], ref_key: &RefKey, params: &'a AlignmentParams) -> std::io::Result<()> {
        for seq in seqs {
            if seq.len() != self.length {
                return Err(std::io::Error::other(format!(
                    "Inconsistent reference lengths for '{reference_id}|{compound_type}'",
                    reference_id = ref_key.reference_id,
                    compound_type = ref_key.compound_type
                )));
            }
        }

        let profiles = seqs.iter().map(|seq| {
            AlignmentProfiles::new(seq, &params.matrix, params.gap_open, params.gap_extend)
                .with_context("Failed to build alignment profiles")
        });
        profiles.process_results(|iter| self.profiles.extend(iter))?;

        Ok(())
    }

    /// Finds the best local Smith-Waterman alignment for a query sequence
    /// against all profiles in this group.
    ///
    /// The alignment with the highest score is considered best, or `None` is
    /// returned if no alignment was found. The `states` in the alignment will
    /// only include `M`, `I`, `D`, and `S`. Furthermore, any alignments are
    /// guaranteed to start and end with `M` states, excluding soft clipping.
    pub fn best_alignment<T: AsRef<[u8]> + ?Sized>(&self, query: &T) -> Option<Alignment<u32>> {
        // TODO: The use of get feels questionable. Even if it is
        // unlikely to ever happen, it feels like it should be an error/warning.
        let mut alignments = self
            .profiles
            .iter()
            .filter_map(|p| p.sw_align_from_i16(query.as_query_src()).get());

        let mut best_alignment = alignments.next()?;

        for alignment in alignments {
            // Instead of using max_by_key, we use manual comparison to ensure
            // the first maximum is returned (not the last)
            if alignment > best_alignment {
                best_alignment = alignment;
            }
        }

        Some(best_alignment)
    }
}
