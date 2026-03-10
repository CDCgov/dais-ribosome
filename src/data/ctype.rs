//! Compound type hierarchy: ctype → reference_id → protein_product.
use super::{
    error::ModuleLoadError,
    exons::Exons,
    keys::{RefKey, SpecKey},
    products::ProductSpec,
    profiles::{AlignmentProfiles, AlignmentWeights, build_alignment_profile},
    refs::ReferenceMap,
    spec::CdsSpecMap,
    weights::CodonWeightMatrix,
};
use std::collections::HashMap;
use zoe::alignment::Alignment;
use zoe::data::nucleotides::Nucleotides;
use zoe::prelude::*;

/// Top-level index: ctype string → compound type data.
pub type CompoundTypeMap<'a> = HashMap<String, Vec<ReferenceGroup<'a>>>;

/// References sharing the same reference_id within a ctype.
#[derive(Debug)]
pub(crate) struct ReferenceGroup<'a> {
    pub(crate) reference_id: String,
    pub(crate) length:       usize,
    pub(crate) profiles:     Vec<AlignmentProfiles<'a>>,
    pub(crate) proteins:     Vec<ProductSpec>,
}

impl<'a> ReferenceGroup<'a> {
    pub fn iter_proteins(&self) -> impl Iterator<Item = &ProductSpec> {
        self.proteins.iter()
    }

    /// Find the best alignment for a query sequence against all profiles in this group.
    ///
    /// Returns the alignment with the highest score, or `None` if no alignment was found.
    pub fn best_alignment<T: AsRef<[u8]> + ?Sized>(&self, query: &T) -> Option<Alignment<u32>> {
        self.profiles
            .iter()
            .filter_map(|p| p.sw_align_from_i16(query.as_query_src()).get())
            .max_by_key(|a| a.score)
    }
}

/// Build the ctype map from raw loaded data.
pub(crate) fn build_ctype_map<'a>(
    references: &'a ReferenceMap, cds_spec: CdsSpecMap, mut codon_weights: CodonWeightMatrix,
    weight_matrices: &'a AlignmentWeights,
) -> Result<CompoundTypeMap<'a>, ModuleLoadError> {
    let mut ctype_refs: HashMap<&str, Vec<(&str, &Vec<Nucleotides>)>> = HashMap::new();
    for (ref_key, seqs) in references {
        ctype_refs
            .entry(&ref_key.compound_type)
            .or_default()
            .push((&ref_key.reference_id, seqs));
    }

    // Group CDS specs by (reference_id, ctype) using RefKey
    let mut spec_by_ref: HashMap<RefKey, Vec<(String, Exons)>> = HashMap::new();
    for (spec_key, ctype_exons) in cds_spec {
        let (ctype, exons) = ctype_exons.into_ctype_exons();
        let key = RefKey::new(spec_key.reference_id, ctype);
        spec_by_ref.entry(key).or_default().push((spec_key.protein_product, exons));
    }

    let mut result = HashMap::new();

    for (ctype, ref_entries) in ctype_refs {
        let (aln_params, weight_matrix) = weight_matrices.get(ctype);

        // Group by reference_id within this ctype
        let mut ref_groups_map: HashMap<&str, Vec<&Nucleotides>> = HashMap::new();
        for (ref_id, seqs) in &ref_entries {
            ref_groups_map.entry(ref_id).or_default().extend(seqs.iter());
        }

        let mut reference_groups = Vec::with_capacity(ref_groups_map.len());

        for (reference_id, seqs) in ref_groups_map {
            // Ensure all sequences are the same length
            let length = seqs.first().map(|s| s.len()).unwrap_or(0);
            for seq in &seqs {
                if seq.len() != length {
                    return Err(ModuleLoadError::validation(format!(
                        "Inconsistent lengths for '{reference_id}|{ctype}'"
                    )));
                }
            }

            // Build the profiles for the equal-length sequences
            let profiles = seqs
                .iter()
                .map(|seq| build_alignment_profile(seq, weight_matrix, aln_params))
                .collect::<Result<Vec<_>, _>>()?;

            let proteins = spec_by_ref
                .remove(&RefKey::new(reference_id, ctype))
                .unwrap_or_default()
                .into_iter()
                .map(|(protein_name, exons)| {
                    let spec_key = SpecKey::new(reference_id, &protein_name);
                    ProductSpec {
                        name: protein_name,
                        exons,
                        codon_weights: codon_weights.remove(&spec_key),
                    }
                })
                .collect();

            reference_groups.push(ReferenceGroup {
                reference_id: reference_id.to_string(),
                length,
                profiles,
                proteins,
            });
        }

        result.insert(ctype.to_string(), reference_groups);
    }

    Ok(result)
}
