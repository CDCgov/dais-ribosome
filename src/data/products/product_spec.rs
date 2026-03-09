use crate::data::{
    exons::Exons,
    keys::CodonKey,
    products::Product,
    ranges::StateRange,
    weights::{CodonPositionWeights, DEFAULT_CODON_STATS},
};

/// Protein product specification.
#[derive(Debug)]
pub(crate) struct ProductSpec {
    /// Protein or peptide name
    pub(crate) name:          String,
    /// Exon Coordinates and translation rules
    pub(crate) exons:         Exons,
    /// The protein product specific codon positive weight matrix
    pub(crate) codon_weights: Option<CodonPositionWeights>,
}

impl ProductSpec {
    pub(crate) fn make_product_ranges<'a>(&'a self, state_ranges: &[StateRange]) -> Product<'a> {
        let mut product_ranges = Vec::with_capacity(state_ranges.len());
        for state in state_ranges.iter() {
            for exon in &self.exons.coords {
                // A state can span multiple exons: e.g., long match for full contig
                // An exon can span multiple states: e.g., a match with an indel
                if let Some(cds_state) = state.intersect_exon(exon) {
                    product_ranges.push(cds_state);
                }
            }
        }

        Product {
            product_ranges,
            product_spec: self,
            stop_extension_query_range: None,
        }
    }

    /// Compare two codons at the same position: returns true if `left >=
    /// right`.
    ///
    /// Uses position-specific weights if available, falling back to
    /// `DEFAULT_CODON_STATS` when both position-specific counts are zero.
    ///
    /// ## Arguments
    ///
    /// - `left` - Left codon (uppercase)
    /// - `right` - Right codon (uppercase)
    /// - `codon_position` - 1-based codon position for weight lookup
    ///
    /// ## Validity
    ///
    /// The `left` and `right` codons must contain unaligned, uppercase IUPAC
    /// bases.
    pub(crate) fn codon_left_ge_right(&self, left: [u8; 3], right: [u8; 3], codon_position: u32) -> bool {
        let (mut x, mut y) = (0u32, 0u32);

        if let Some(w) = &self.codon_weights {
            x = w.get(&CodonKey::new(codon_position, left)).copied().unwrap_or(0);
            y = w.get(&CodonKey::new(codon_position, right)).copied().unwrap_or(0);
        }

        // Fall back to default stats if both position-specific counts are zero
        if x == 0 && y == 0 {
            x = DEFAULT_CODON_STATS.get(&left).copied().unwrap_or(0);
            y = DEFAULT_CODON_STATS.get(&right).copied().unwrap_or(0);
        }

        x >= y
    }

    /// Compares the same codon at two different positions: returns true if
    /// `pos_left >= pos_right`.
    ///
    /// Used for deletion frame correction where we need to decide which
    /// position the pivot codon should be assigned to.
    ///
    /// ## Arguments
    ///
    /// - `pos_left` - 1-based left codon position
    /// - `pos_right` - 1-based right codon position
    /// - `codon` - The codon to compare (uppercase)
    /// - `preference` - Which direction to prefer when both counts are zero
    pub(crate) fn codon_pos_left_ge_right(
        &self, pos_left: u32, pos_right: u32, codon: [u8; 3], preference: ShiftPreference,
    ) -> bool {
        let (mut x, mut y) = (0u32, 0u32);

        if let Some(w) = &self.codon_weights {
            x = w.get(&CodonKey::new(pos_left, codon)).copied().unwrap_or(0);
            y = w.get(&CodonKey::new(pos_right, codon)).copied().unwrap_or(0);
        }

        // Use preference to break ties when both are zero
        if x == 0 && y == 0 {
            return preference == ShiftPreference::Left;
        }

        x >= y
    }
}

/// Direction preference for tie-breaking when comparing codon positions.
#[derive(Clone, Copy, PartialEq, Eq)]
pub(crate) enum ShiftPreference {
    /// Prefer left shift on tie
    Left,
    /// Prefer right shift on tie
    Right,
}
