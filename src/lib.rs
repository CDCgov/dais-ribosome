mod annotation;
pub(crate) mod config;
pub(crate) mod data;
pub mod errors;
pub(crate) mod hashing;
pub mod outputs;
pub mod tsv;

pub use config::{annotation_module::AnnotationModule, toml};
pub use data::{NoNucleotides, QueryRecord, ranges};

use std::ops::ControlFlow;
use zoe::{alignment::AlignmentStates, data::cigar::Ciglet};

trait IteratorExt: Iterator {
    fn take_until_inclusive<F>(self, f: F) -> TakeUntilInclusive<Self, F>
    where
        Self: Sized,
        F: FnMut(&Self::Item) -> bool, {
        TakeUntilInclusive::new(self, f)
    }
}

impl<I> IteratorExt for I where I: Iterator {}

/// An iterator adaptor that consumes elements until the given predicate is
/// `true`, including that element.
///
/// This is based on `TakeWhileInclusive` from Itertools, but with the predicate
/// negated.
#[derive(Clone)]
struct TakeUntilInclusive<I, F> {
    iter:      I,
    predicate: F,
    done:      bool,
}

impl<I, F> TakeUntilInclusive<I, F>
where
    I: Iterator,
    F: FnMut(&I::Item) -> bool,
{
    /// Create a new [`TakeUntilInclusive`] from an iterator and a predicate.
    pub(crate) fn new(iter: I, predicate: F) -> Self {
        Self {
            iter,
            predicate,
            done: false,
        }
    }
}

impl<I, F> Iterator for TakeUntilInclusive<I, F>
where
    I: Iterator,
    F: FnMut(&I::Item) -> bool,
{
    type Item = I::Item;

    fn next(&mut self) -> Option<Self::Item> {
        if self.done {
            None
        } else {
            self.iter.next().inspect(|item| {
                if (self.predicate)(item) {
                    self.done = true;
                }
            })
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        if self.done {
            (0, Some(0))
        } else {
            (0, self.iter.size_hint().1)
        }
    }

    fn fold<B, Fold>(mut self, init: B, mut f: Fold) -> B
    where
        Fold: FnMut(B, Self::Item) -> B, {
        if self.done {
            init
        } else {
            let out = self.iter.try_fold(init, |mut acc, item| {
                let exit = (self.predicate)(&item);
                acc = f(acc, item);
                if exit {
                    ControlFlow::Break(acc)
                } else {
                    ControlFlow::Continue(acc)
                }
            });

            match out {
                ControlFlow::Continue(acc) | ControlFlow::Break(acc) => acc,
            }
        }
    }
}

/// An extension trait for [`AlignmentStates`], providing functionality custom
/// to DAIS-ribosome.
trait AlignmentStatesExt {
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
