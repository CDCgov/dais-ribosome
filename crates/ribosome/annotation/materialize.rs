//! Procedures for converting range-based [`Product`] outputs into materialized
//! [`ComputedProduct`].

use crate::{
    IteratorExt, QueryRecord,
    data::ranges::{CdsDeletionRange, CdsInsertionRange, CdsMatchRange, CdsStateRange},
    hashing::{nt_id_iupac, variant_hash_iupac},
    outputs::{ComputedDeletion, ComputedInsertion, ComputedProduct, Product},
    ranges::CdsCoord,
};
use std::ops::Range;
use zoe::prelude::{AminoAcids, Len, Nucleotides, Slice, Translate};

impl<'a> Product<'a> {
    /// Computes the output data for this product, materializing all ranges into
    /// sequences using `query`.
    ///
    /// The following cases describe the length of the output sequences:
    ///
    /// - If `product.product_ranges` is empty, then the sequence and coordinate
    ///   fields will be empty. `trailing_cds_unaligned` will be non-zero.
    /// - If `product.product_ranges` contains only a deletion, then `cds_aln`
    ///   and `aa_aln` will be non-empty. `cds_seq`, `aa_seq`, and the
    ///   coordinate fields will be empty.
    /// - Otherwise, all the sequence and coordinate fields will be non-empty.
    ///
    /// ## Validity
    ///
    /// The `query` must be the same record used to form the product.
    pub fn materialize(&self, query: &QueryRecord) -> ComputedProduct<'a> {
        // Compute all the fields that rely on incremental updates until the
        // first stop codon
        let incremental = ComputedIncrementalProducts::new(query, self);

        let ComputedIncrementalProducts {
            cds_aln,
            aa_aln,
            cds_seq,
            has_insertion,
            has_shift_indel,
            query_coords,
            cds_coords,
            insertions,
            deletions,
            trailing_cds_unaligned,
        } = incremental;

        // Form aa_seq by splicing insertions into aa_aln (and removing
        // deletions)
        let aa_seq = {
            let mut out = AminoAcids::new();

            let mut aa_aln_iter = aa_aln.iter().copied();

            // The number of amino acids consumed from aa_aln so far
            let mut num_consumed = 0;
            for insertion in &insertions {
                // 1-based index after which insertion occurs is equivalently
                // the count of the number of amino acids before the insertion.
                let num_to_consume = insertion.upstream_aa_pos - num_consumed;

                // Extend with non-deleted aligned residues before the insertion
                let aligned_before_insertion = aa_aln_iter.by_ref().take(num_to_consume).filter(|&b| b != b'-' && b != b'.');

                out.extend(aligned_before_insertion);

                // Validity: ComputedInsertion::inserted_aa is uppercase IUPAC
                out.extend_from_slice(&insertion.inserted_aa);

                num_consumed += num_to_consume;
            }

            // Extend with non-deleted aligned residues after the last insertion
            out.extend(aa_aln_iter.filter(|&b| b != b'-' && b != b'.'));

            out
        };

        // Validity: cds_seq contains unaligned uppercase IUPAC
        let cds_id = nt_id_iupac(&cds_seq);

        // Validity: aa_seq contains unaligned uppercase IUPAC (since it is
        // derived from StdGeneticCode with gaps filtered)
        let variant_hash = variant_hash_iupac(&aa_seq);

        ComputedProduct {
            name: &self.product_spec.name,
            cds_seq,
            cds_aln,
            cds_id,
            aa_seq,
            aa_aln,
            variant_hash,
            has_insertion,
            has_shift_indel,
            query_coords,
            cds_coords,
            insertions,
            deletions,
            trailing_cds_unaligned,
        }
    }
}

/// All fields of [`ComputedProduct`] that are populated incrementally via
/// iterating over the `product_ranges`.
///
/// None of the fields are populated with any ranges or residues past the first
/// encountered stop codon. From these fields, the rest of [`ComputedProduct`]
/// can be computed.
///
/// [`ComputedProduct`]: crate::outputs::ComputedProduct
struct ComputedIncrementalProducts {
    /// See [`ComputedProduct::cds_aln`]. This is populated from
    /// [`IncrementalAccumulator::cds_aln`].
    cds_aln:                Nucleotides,
    /// See [`ComputedProduct::aa_aln`]. This is populated from
    /// [`IncrementalAccumulator::aa_aln`].
    aa_aln:                 AminoAcids,
    /// See [`ComputedProduct::cds_seq`]. This is populated from
    /// [`DependentFields::cds_seq`].
    cds_seq:                Nucleotides,
    /// See [`ComputedProduct::has_insertion`]. This is populated from
    /// [`DependentFields::has_insertion`].
    has_insertion:          bool,
    /// See [`ComputedProduct::has_shift_indel`]. This is populated from
    /// [`DependentFields::has_insertion`].
    has_shift_indel:        bool,
    /// See [`ComputedProduct::query_coords`]. This is populated from
    /// [`DependentFields::query_coords`].
    query_coords:           Vec<Range<usize>>,
    /// See [`ComputedProduct::cds_coords`]. This is populated from
    /// [`DependentFields::cds_coords`].
    cds_coords:             Vec<CdsCoord>,
    /// See [`ComputedProduct::insertions`]. This is populated from
    /// [`DependentFields::insertions`].
    insertions:             Vec<ComputedInsertion>,
    /// See [`ComputedProduct::deletions`]. This is populated from
    /// [`DependentFields::deletions`].
    deletions:              Vec<ComputedDeletion>,
    /// See [`ComputedProduct::trailing_cds_unaligned`].
    trailing_cds_unaligned: usize,
}

impl ComputedIncrementalProducts {
    /// Computes the incremental products from a `query` and `product`.
    ///
    /// The following cases describe the length of the output sequences:
    ///
    /// - If `product.product_ranges` is empty, then the sequence and coordinate
    ///   fields will be empty. `trailing_cds_unaligned` will be non-zero.
    /// - If `product.product_ranges` contains only a deletion, then `cds_aln`
    ///   and `aa_aln` will be non-empty, while `cds_seq` and the coordinate
    ///   fields will be empty.
    /// - Otherwise, all the sequence and coordinate fields will be non-empty.
    fn new(query: &QueryRecord, product: &Product) -> Self {
        let range_capacity = product.product_ranges.len();

        let mut out = IncrementalAccumulator::new(product.leading_cds_unaligned, range_capacity);

        let end_cds_index = out.populate_from(query, product);

        let trailing_cds_unaligned = product.product_spec.exons.cds_len() - end_cds_index;

        Self {
            cds_aln: out.cds_aln,
            aa_aln: out.aa_aln,
            cds_seq: out.dependent_fields.cds_seq,
            has_insertion: out.dependent_fields.has_insertion,
            has_shift_indel: out.dependent_fields.has_shift_indel,
            query_coords: out.dependent_fields.query_coords,
            cds_coords: out.dependent_fields.cds_coords,
            insertions: out.dependent_fields.insertions,
            deletions: out.dependent_fields.deletions,
            trailing_cds_unaligned,
        }
    }
}

/// A private struct for aiding in the construction of
/// [`ComputedIncrementalProducts`] which supports incremental updates for each
/// range.
struct IncrementalAccumulator {
    /// The field for incrementally building [`ComputedProduct::cds_aln`].
    ///
    /// This is updated eagerly on each call to [`extend_from_match`] or
    /// [`extend_from_deletion`].
    ///
    /// ## Validity
    ///
    /// This must only contain uppercase IUPAC, padding `.`, and gaps `-`. Both
    /// `U` and `T` are allowed.
    ///
    /// [`extend_from_match`]: IncrementalAccumulator::extend_from_match
    /// [`extend_from_deletion`]: IncrementalAccumulator::extend_from_deletion
    cds_aln: Nucleotides,

    /// The field for incrementally building [`ComputedProduct::aa_aln`].
    ///
    /// This is updated from `cds_aln` on each call to [`extend_from_match`] or
    /// when accumulation is complete. This is built incrementally to facilitate
    /// detection of the stop codon.
    ///
    /// ## Validity
    ///
    /// This must only contain uppercase IUPAC, partial codons, stop codons,
    /// padding `.`, and gaps `-`.
    ///
    /// [`extend_from_match`]: IncrementalAccumulator::extend_from_match
    aa_aln: AminoAcids,

    /// The index in `cds_aln` marking the start of the bases which have not yet
    /// been translated into `aa_aln`.
    ///
    /// If `untranslated_start == cds_aln.len()`, then all bases have been
    /// translated so far.
    untranslated_start: usize,

    /// All computed product fields which are dependent on the aligned
    /// sequences, and hence will be truncated if they are truncated.
    dependent_fields: DependentFields,
}

impl IncrementalAccumulator {
    /// Initializes an [`IncrementalAccumulator`] containing just the gap
    /// indicated by `leading_cds_unaligned`.
    fn new(leading_cds_unaligned: usize, range_capacity: usize) -> Self {
        Self {
            cds_aln:            Nucleotides::from(vec![b'.'; leading_cds_unaligned]),
            aa_aln:             AminoAcids::new(),
            untranslated_start: 0,
            dependent_fields:   DependentFields::new(range_capacity),
        }
    }

    /// Populates the [`IncrementalAccumulator`] from all the ranges in
    /// `product`, returning the exclusive-end index of the accumulated
    /// sequences within the coding sequence (used in initializing
    /// `trailing_cds_unaligned`).
    ///
    /// The following cases describe the length of the output sequences:
    ///
    /// - If `product.product_ranges` is empty, then the sequence and coordinate
    ///   fields will be empty.
    /// - If `product.product_ranges` contains only a deletion, then `cds_aln`
    ///   and `aa_aln` will be non-empty, while `cds_seq` and the coordinate
    ///   fields will be empty.
    /// - Otherwise, all the sequence and coordinate fields will be non-empty.
    fn populate_from(&mut self, query: &QueryRecord, product: &Product) -> usize {
        let ranges = product.product_ranges.iter().enumerate();

        // Populate from ranges up to the last
        for (i, state) in ranges {
            // Validity: product_ranges.len() - 1 will not underflow since it is
            // non-empty (otherwise the loop wouldn't be running)
            let is_terminal = i == 0 || i == product.product_ranges.len() - 1;

            match state {
                CdsStateRange::M(m) => {
                    // Short circuit if a stop codon is found
                    if let Some(end_cds_index) = self.extend_from_match(query, m) {
                        return end_cds_index;
                    }
                }
                CdsStateRange::I(ins) => self.extend_from_insertion(query, ins),
                CdsStateRange::D(del) => self.extend_from_deletion(del, is_terminal),
            };
        }

        // No stop codon reached, so finish the translation
        self.update_translation();

        // Possibly add a partial codon
        if self.untranslated_start < self.cds_aln.len() {
            self.aa_aln.push(b'~');
        }

        product
            .product_ranges
            .iter()
            .rev()
            .find_map(|s| match s {
                CdsStateRange::M(m) => Some(m.cds_range.end),
                CdsStateRange::D(d) => Some(d.cds_range.end),
                _ => None,
            })
            .unwrap_or(0)
    }

    /// Updates the accumulator with a match range. If a stop codon is reached,
    /// translation halts and the end index in the coding sequence is returned.
    ///
    /// This method guarantees that `cds_aln`, `cds_seq`, `query_coords`, and
    /// `cds_coords` will grow in size.
    ///
    /// ## Validity
    ///
    /// `range` must have a non-zero length.
    #[must_use]
    fn extend_from_match(&mut self, query: &QueryRecord, range: &CdsMatchRange) -> Option<usize> {
        let slice = &query.nucleotides()[range.query_range.clone()];
        // Validity: QueryRecord only contains uppercase IUPAC
        self.cds_aln.extend_from_slice(slice);

        self.update_translation();

        // Check if we reached the stop codon
        if self.aa_aln.ends_with(b"*") {
            let range_until_stop = {
                let untranslated_len = self.cds_aln.len() - self.untranslated_start;
                let mut r = range.clone();
                r.cut_end(untranslated_len);
                r
            };

            // Shrink the aligned coding sequence
            self.cds_aln.shorten_to(self.untranslated_start);

            self.dependent_fields.extend_from_match(query, &range_until_stop);

            Some(range_until_stop.cds_range.end)
        } else {
            self.dependent_fields.extend_from_match(query, range);
            None
        }
    }

    /// Updates the accumulator with an insertion range.
    ///
    /// This method guarantees that `cds_seq`, `query_coords`, and `cds_coords`
    /// will grow in size.
    ///
    /// ## Validity
    ///
    /// `range` must have a non-zero length.
    fn extend_from_insertion(&mut self, query: &QueryRecord, range: &CdsInsertionRange) {
        self.dependent_fields.extend_from_insertion(query, range);
    }

    /// Updates the accumulator with a deletion range.
    ///
    /// This method guarantees that `cds_aln` will grow in size.
    ///
    /// ## Validity
    ///
    /// `range` must have a non-zero length.
    fn extend_from_deletion(&mut self, range: &CdsDeletionRange, is_terminal: bool) {
        self.cds_aln.extend(std::iter::repeat_n(b'-', range.len()));
        self.dependent_fields.extend_from_deletion(range, is_terminal);
    }

    /// Updates the `aa_aln` translation to include as many of the bytes in
    /// `cds_aln` as possible. Only complete (length-3) codons are translated.
    fn update_translation(&mut self) {
        let untranslated = self.cds_aln.slice(self.untranslated_start..);
        self.aa_aln
            .extend(untranslated.to_aa_iter_exact().take_until_inclusive(|aa| *aa == b'*'));
        self.untranslated_start = self.aa_aln.len() * 3;
    }
}

/// All incrementally computed product fields which are dependent on the aligned
/// sequences, and hence will be truncated if they are truncated.
struct DependentFields {
    /// The field for incrementally building [`ComputedProduct::cds_seq`].
    ///
    /// ## Validity
    ///
    /// This must only contain unaligned uppercase IUPAC. Both `U` and `T` are
    /// allowed.
    cds_seq: Nucleotides,

    /// The boolean flag for [`ComputedProduct::has_insertion`], which starts as
    /// `false` and gets set to `true` when an insertion is encountered.
    has_insertion: bool,

    /// The boolean flag for [`ComputedProduct::has_shift_indel`], which starts
    /// as `false` and gets set to `true` when a shift indel is encountered.
    has_shift_indel: bool,

    /// The field for incrementally building [`ComputedProduct::query_coords`].
    ///
    /// ## Validity
    ///
    /// This must be the same length as `cds_coords`.
    query_coords: Vec<Range<usize>>,

    /// The field for incrementally building [`ComputedProduct::cds_coords`].
    ///
    /// ## Validity
    ///
    /// This must be the same length as `query_coords`.
    cds_coords: Vec<CdsCoord>,

    /// The field for incrementally collecting [`ComputedProduct::insertions`].
    insertions: Vec<ComputedInsertion>,

    /// The field for incrementally collecting [`ComputedProduct::deletions`].
    deletions: Vec<ComputedDeletion>,
}

impl DependentFields {
    /// Initializes an empty [`DependentFields`] with the given number of ranges
    /// as a starting capacity for `query_coords` and `cds_coords` fields.
    fn new(range_capacity: usize) -> Self {
        Self {
            cds_seq:         Nucleotides::new(),
            has_insertion:   false,
            has_shift_indel: false,
            query_coords:    Vec::with_capacity(range_capacity),
            cds_coords:      Vec::with_capacity(range_capacity),
            insertions:      Vec::new(),
            deletions:       Vec::new(),
        }
    }

    /// Extends the dependent fields from a possibly-truncated match range.
    ///
    /// This method guarantees that `cds_seq`, `query_coords`, and `cds_coords`
    /// will grow in size.
    ///
    /// ## Validity
    ///
    /// `range` must have a non-zero length.
    fn extend_from_match(&mut self, query: &QueryRecord, range: &CdsMatchRange) {
        let slice = &query.nucleotides()[range.query_range.clone()];

        // Validity: slice is from QueryRecord, hence meets the requirements
        self.cds_seq.extend_from_slice(slice);

        // Validity: we push to query_coords and cds_coords at the same time
        self.query_coords.push(range.query_range.clone());
        self.cds_coords.push(CdsCoord::M(range.cds_range.clone()));
    }

    /// Extends the dependent fields from a insertion range.
    ///
    /// This method guarantees that `cds_seq`, `query_coords`, and `cds_coords`
    /// will grow in size.
    ///
    /// ## Validity
    ///
    /// `range` must have a non-zero length.
    fn extend_from_insertion(&mut self, query: &QueryRecord, range: &CdsInsertionRange) {
        let slice = &query.nucleotides()[range.query_range.clone()];

        // Validity: slice is from QueryRecord
        let (computed_insertion, filtered) = ComputedInsertion::new(range.cds_index, slice);

        // Validity: ComputedInsertion::inserted_nt meets validity requirements
        self.cds_seq.extend_from_slice(&computed_insertion.inserted_nt);

        if !filtered {
            self.insertions.push(computed_insertion);
        }
        self.has_insertion = true;

        // Validity: we push to query_coords and cds_coords at the same time
        self.query_coords.push(range.query_range.clone());
        self.cds_coords.push(CdsCoord::I(range.cds_index));

        if !range.len().is_multiple_of(3) {
            self.has_shift_indel = true;
        }
    }

    /// Extends the dependent fields from a deletion range.
    ///
    /// ## Validity
    ///
    /// `range` must have a non-zero length.
    fn extend_from_deletion(&mut self, range: &CdsDeletionRange, is_terminal: bool) {
        // Only allow has_shift_indel to be updated for a deletion if it is not
        // at the beginning or end of the ranges for the given exon
        if !is_terminal && !range.len().is_multiple_of(3) {
            self.has_shift_indel = true;
        }

        let in_frame = range.cds_range.start.is_multiple_of(3) && range.len().is_multiple_of(3);

        // The 1-based inclusive start. Floor divide so that
        // deleting any part of a codon deletes the amino acid. Add
        // 1 to make it 1-based.
        let del_aa_start = (range.cds_range.start / 3) + 1;

        // The 1-based inclusive end. Ceiling divide so that
        // deleting any part of a codon deletes the amino acid.
        // 0-based exclusive end is equivalent to 1-based inclusive
        // end.
        let del_aa_end = range.cds_range.end.div_ceil(3);

        let del_cds_len = range.len();
        let del_aa_len = del_cds_len.div_ceil(3);

        self.deletions.push(ComputedDeletion {
            del_aa_start,
            del_aa_end,
            del_aa_len,
            in_frame,
            del_cds_start: range.cds_range.start + 1,
            del_cds_end: range.cds_range.end,
            del_cds_len,
        });
    }
}
