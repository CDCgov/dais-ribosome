//! Compound type hierarchy: ctype → reference_id → protein_product.
use super::{
    error::ModuleLoadError,
    keys::{RefKey, SpecKey},
    products::ProductSpec,
    refs::ReferenceMap,
    weights::CodonWeightMatrix,
};
use crate::{
    config::toml::{AlignmentParams, AlignmentWeights},
    data::{exons::Exons, spec::CdsSpecMap},
};
use std::collections::HashMap;
use zoe::alignment::{Alignment, SharedProfiles};
use zoe::prelude::*;
use zoe::{data::nucleotides::Nucleotides, iter_utils::ProcessResultsExt};

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
    pub(crate) proteins:     Vec<ProductSpec>,
}

impl<'a> ReferenceGroup<'a> {
    pub fn new(
        ref_key: &RefKey, seqs: &'a [Nucleotides], params: &'a AlignmentParams,
        cds_spec: &mut HashMap<RefKey, Vec<(String, Exons)>>, codon_weights: &mut CodonWeightMatrix,
    ) -> Result<Self, ModuleLoadError> {
        let length = seqs.first().map_or(0, Nucleotides::len);
        for seq in seqs {
            if seq.len() != length {
                return Err(ModuleLoadError::validation(format!(
                    "Inconsistent lengths for '{reference_id}|{compound_type}'",
                    reference_id = ref_key.reference_id,
                    compound_type = ref_key.compound_type
                )));
            }
        }

        let profiles = seqs
            .iter()
            .map(|seq| AlignmentProfiles::new(seq, &params.matrix, params.gap_open, params.gap_extend))
            .collect::<Result<Vec<_>, _>>()?;

        let proteins = cds_spec
            .remove(ref_key)
            .unwrap_or_default()
            .into_iter()
            .map(|(protein_name, exons)| {
                let spec_key = SpecKey::new(&ref_key.reference_id, &protein_name);
                ProductSpec {
                    name: protein_name,
                    exons,
                    codon_weights: codon_weights.remove(&spec_key),
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

    pub fn extend(
        &mut self, seqs: &'a [Nucleotides], ref_key: &RefKey, params: &'a AlignmentParams,
    ) -> Result<(), ModuleLoadError> {
        for seq in seqs {
            if seq.len() != self.length {
                return Err(ModuleLoadError::validation(format!(
                    "Inconsistent lengths for '{reference_id}|{compound_type}'",
                    reference_id = ref_key.reference_id,
                    compound_type = ref_key.compound_type
                )));
            }
        }

        let profiles = seqs
            .iter()
            .map(|seq| AlignmentProfiles::new(seq, &params.matrix, params.gap_open, params.gap_extend));
        profiles.process_results(|iter| self.profiles.extend(iter))?;

        Ok(())
    }

    /// Finds the best alignment for a query sequence against all profiles in
    /// this group.
    ///
    /// Returns the alignment with the highest score, or `None` if no alignment
    /// was found.
    pub fn best_alignment<T: AsRef<[u8]> + ?Sized>(&self, query: &T) -> Option<Alignment<u32>> {
        // TODO: Why did we hard code from i16? Will this potentially disagree
        // with aligner? Also, the use of get feels questionable. Even if it is
        // unlikely to ever happen, it feels like it should be an error/warning.
        self.profiles
            .iter()
            .filter_map(|p| p.sw_align_from_i16(query.as_query_src()).get())
            .max_by_key(|a| a.score)
    }
}

/// Build the ctype map from raw loaded data.
pub(crate) fn build_ctype_map<'a>(
    references: &'a ReferenceMap, mut cds_spec: CdsSpecMap, mut codon_weights: CodonWeightMatrix,
    alignment_weights: &'a AlignmentWeights,
) -> Result<CompoundTypeMap<'a>, ModuleLoadError> {
    let mut ctype_map: HashMap<String, Vec<ReferenceGroup<'a>>> = HashMap::new();

    for (ref_key, seqs) in references {
        let params = alignment_weights.get(&ref_key.compound_type);

        // Get the list of groups for the given compound type
        let groups = ctype_map.entry(ref_key.compound_type.to_string()).or_default();

        // See if there is an existing entry in the list of groups for the given
        // reference ID
        if let Some(group) = groups.iter_mut().find(|group| group.reference_id == ref_key.reference_id) {
            // Update that reference group
            group.extend(seqs, ref_key, params)?;
        } else {
            // Add a new reference group
            groups.push(ReferenceGroup::new(ref_key, seqs, params, &mut cds_spec, &mut codon_weights)?);
        }
    }

    Ok(ctype_map)
}
