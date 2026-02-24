#![feature(int_format_into)]

// =============================================================================
// Module hierarchy
// =============================================================================

use zoe::{alignment::AlignmentStates, data::cigar::Ciglet};

/// Configuration loading and path resolution.
pub mod config;

/// Data loading and structures for module resources.
pub mod data;

/// Protein annotation engine.
pub mod annotation;

/// An extension trait for [`AlignmentStates`], providing functionality custom
/// to DAIS-ribosome.
pub trait AlignmentStatesExt {
    /// Extends an [`AlignmentStates`] to the left by prepending `inc` match
    /// states, replacing corresponding increments of soft clipping if present.
    fn extend_left(&mut self, inc: usize);

    /// Extends an [`AlignmentStates`] to the right by appending `inc` match
    /// states, replacing corresponding increments of soft clipping if present.
    fn extend_right(&mut self, inc: usize);
}

impl AlignmentStatesExt for AlignmentStates {
    fn extend_left(&mut self, inc: usize) {
        // Is the first ciglet soft clipping, in which case we have to remove
        // some of it? If not, just prepend as normal.
        if let Some((first, rest)) = self.as_mut_slice().split_first_mut()
            && first.op == b'S'
        {
            let previous_soft_clipping = first.inc;
            let new_soft_clipping = previous_soft_clipping.saturating_sub(inc);

            // Is the ciglet after the soft clipping a match state? If so, we
            // need to merge into it. If not, we can replace the first ciglet
            // with our new one, and then optionally add back soft clipping.
            if let Some(second) = rest.first_mut()
                && second.op == b'M'
            {
                second.inc += inc;
                // Is the soft clipping remaining positive? If so, mutate the
                // first ciglet. If not, we need to remove it entirely.
                if new_soft_clipping > 0 {
                    first.inc = new_soft_clipping;
                } else {
                    self.as_mut_vec().remove(0);
                }
            } else {
                *first = Ciglet { inc, op: b'M' };
                self.prepend_soft_clip(new_soft_clipping);
            }
        } else {
            self.prepend_inc_op(inc, b'M');
        }
    }

    fn extend_right(&mut self, inc: usize) {
        // Is the last ciglet soft clipping, in which case we have to remove
        // some of it? If not, just prepend as normal.
        if let Some((last, rest)) = self.as_mut_slice().split_last_mut()
            && last.op == b'S'
        {
            let previous_soft_clipping = last.inc;
            let new_soft_clipping = previous_soft_clipping.saturating_sub(inc);

            // Is the ciglet before the soft clipping a match state? If so, we
            // need to merge into it. If not, we can replace the last ciglet
            // with our new one, and then optionally add back soft clipping.
            if let Some(second_to_last) = rest.last_mut()
                && second_to_last.op == b'M'
            {
                second_to_last.inc += inc;
                // Is the soft clipping remaining positive? If so, mutate the
                // last ciglet. If not, we need to remove it entirely.
                if new_soft_clipping > 0 {
                    last.inc = new_soft_clipping;
                } else {
                    self.as_mut_vec().pop();
                }
            } else {
                *last = Ciglet { inc, op: b'M' };
                self.soft_clip(new_soft_clipping);
            }
        } else {
            self.add_inc_op(inc, b'M');
        }
    }
}
