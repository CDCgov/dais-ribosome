use zoe::prelude::Len;

use crate::{
    config::toml::Formatting,
    data::{
        Nullable,
        products::{ComputedDeletion, ComputedInsertion, ComputedProduct, Product},
    },
};
use std::fmt::{self, Display};

// Can remove later with the regression testing feature later
#[allow(dead_code)]
/// A single row for `.seq` output.
pub struct SeqRow<'a> {
    pub id: &'a str,
    pub ctype: &'a str,
    pub ref_id: &'a str,
    pub(crate) product: &'a Product<'a>,
    pub(crate) computed_product: &'a ComputedProduct,
    pub(crate) formatting: &'a Formatting,
}

impl Display for SeqRow<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let d = self.computed_product;

        let trailing = self.product.trailing_cds_gap_len();

        // Regression: always pad CDS, never pad AA, synthesize empty AA alignments.
        #[cfg(feature = "regression-testing")]
        let pad = ".".repeat(trailing);
        #[cfg(feature = "regression-testing")]
        let cds_rpad = pad.as_str();
        #[cfg(feature = "regression-testing")]
        let aa_rpad = "";
        #[cfg(feature = "regression-testing")]
        let aa_aln = if d.aa_aln.is_empty() {
            vec![b'.'; (d.cds_aln.len() + trailing) / 3].into()
        } else {
            d.aa_aln.clone()
        };

        // Normal: empty data is strictly nullable; skip padding when there is no data.
        #[cfg(not(feature = "regression-testing"))]
        let pad = if d.cds_aln.is_empty() {
            String::new()
        } else {
            ".".repeat(trailing)
        };
        #[cfg(not(feature = "regression-testing"))]
        let cds_rpad = if self.formatting.right_pad_cds && !pad.is_empty() {
            &pad[..trailing]
        } else {
            ""
        };
        #[cfg(not(feature = "regression-testing"))]
        let aa_rpad = if self.formatting.right_pad_aa && !d.aa_aln.is_empty() {
            &pad[..trailing / 3]
        } else {
            ""
        };
        #[cfg(not(feature = "regression-testing"))]
        let aa_aln = Nullable(&d.aa_aln);

        write!(
            f,
            concat!(
                "{id}\t{ctype}\t{ref_id}\t{prot}\t{vh}",
                "\t{aa_seq}\t{aa_aln}{aa_rpad}",
                "\t{cds_id}\t{ins}\t{shift}",
                "\t{cds_seq}\t{cds_aln}{c_rpad}",
                "\t{q_coords}\t{c_coords}"
            ),
            id = self.id,
            ctype = self.ctype,
            ref_id = self.ref_id,
            prot = self.product.product_spec.name,
            vh = Nullable(&d.variant_hash),
            aa_seq = Nullable(&d.aa_seq),
            aa_aln = aa_aln,
            aa_rpad = aa_rpad,
            cds_id = Nullable(&d.cds_id),
            ins = d.has_insertion,
            shift = d.has_shift_indel,
            cds_seq = Nullable(&d.cds_seq),
            cds_aln = d.cds_aln,
            c_rpad = cds_rpad,
            q_coords = Nullable(&d.query_coords),
            c_coords = Nullable(&d.cds_coords),
        )
    }
}

/// A single row for `.ins` output.
pub struct InsRow<'a> {
    pub id:               &'a str,
    pub ctype:            &'a str,
    pub ref_id:           &'a str,
    pub protein:          &'a str,
    pub(crate) insertion: &'a ComputedInsertion,
}

impl Display for InsRow<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let ins = self.insertion;

        #[cfg(feature = "regression-testing")]
        let nts = ins.inserted_nucleotides.to_string().to_lowercase();

        #[cfg(not(feature = "regression-testing"))]
        let nts = &ins.inserted_nucleotides;

        write!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.id,
            self.ctype,
            self.ref_id,
            self.protein,
            ins.upstream_aa,
            nts,
            ins.inserted_residues,
            ins.upstream_nt,
            ins.codon_shift,
        )
    }
}

/// A single row for `.del` output.
pub struct DelRow<'a> {
    pub id: &'a str,
    pub ctype: &'a str,
    pub ref_id: &'a str,
    pub protein: &'a str,
    pub(crate) computed_product: &'a ComputedProduct,
    pub(crate) deletion: &'a ComputedDeletion,
}

impl Display for DelRow<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let del = self.deletion;
        let d = self.computed_product;

        write!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.id,
            self.ctype,
            self.ref_id,
            self.protein,
            Nullable(&d.variant_hash),
            del.del_aa_start,
            del.del_aa_end,
            del.del_aa_len,
            del.in_frame,
            Nullable(&d.cds_id),
            del.del_cds_start,
            del.del_cds_end,
            del.del_cds_len,
        )
    }
}
