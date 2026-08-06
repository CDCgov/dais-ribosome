/// The state at a given index, along with flanking states if they are in
/// bounds.
///
/// This is returned by [`get_state_with_flanking`], and used by [`fix_frame`].
///
/// ## Parameters
///
/// - `'a`: The lifetime of the mutable reference
/// - `T`: The type of the range, such as [`CdsStateRange`] or [`StateRange`]
pub struct StateWithFlanking<'a, T> {
    /// The state two to the left of the current one (`idx-2`)
    pub left2:   Option<&'a mut T>,
    /// The state left of the current one (`idx-1`)
    pub left1:   Option<&'a mut T>,
    /// The current states (`idx`)
    pub current: &'a mut T,
    /// The state right of the current one (`idx+1`)
    pub right1:  Option<&'a mut T>,
    /// The state two to the right of the current one (`idx+2`)
    pub right2:  Option<&'a mut T>,
}

/// Gets the state at `idx` within `ranges`, as well as the two flanking states
/// if available.
///
/// This is a helper function for [`fix_frames`]. If `None` is returned, then
/// the index is out of bounds.
///
/// [`fix_frames`]: Product::fix_frames
#[must_use]
pub fn get_state_with_flanking<T>(idx: usize, product_ranges: &mut [T]) -> Option<StateWithFlanking<'_, T>> {
    let (left, current_and_right) = product_ranges.split_at_mut_checked(idx)?;
    let (current, right) = current_and_right.split_first_mut()?;

    let (left2, left1) = match left {
        [.., left2, left1] => (Some(left2), Some(left1)),
        [left1] => (None, Some(left1)),
        [] => (None, None),
    };

    let (right1, right2) = match right {
        [right1, right2, ..] => (Some(right1), Some(right2)),
        [right1] => (Some(right1), None),
        [] => (None, None),
    };

    Some(StateWithFlanking {
        left2,
        left1,
        current,
        right1,
        right2,
    })
}

pub fn rewrite<T, F>(ranges: &mut Vec<T>, mut f: F)
where
    F: FnMut(StateWithFlanking<T>) -> IdxAdjustment<T>, {
    // The index of the current CdsStateRange to correct/handle
    let mut idx = 0;

    while let Some(states) = get_state_with_flanking(idx, ranges) {
        // Perform any frame fixing on range, which then returns whether to
        // advance the index or not, as well as any states to remove.
        let IdxAdjustment { advance, edits } = f(states);

        if let Some(edits) = edits {
            // Remove the specified states, returning the resulting shift
            // that will need to be applied to idx
            let idx_shift = edits.apply(ranges, idx);

            // Advance idx first, to avoid underflow
            if advance {
                idx += 1;
            }

            // Apply the shift to the index, using wrapping for efficiency
            // (overflow is extremely unlikely)
            idx = idx.wrapping_add_signed(idx_shift);
        } else if advance {
            // Advance idx
            idx += 1;
        }
    }
}

/// The return value for a frame-fixing step (the closure in [`rewrite`]),
/// indicating which states need to be removed, and whether or not the index
/// should be advanced in [`rewrite`].
pub struct IdxAdjustment<T> {
    /// Whether to advance the index, or whether to rehandle the same state
    pub advance: bool,
    /// Any states to remove or insert
    pub edits:   Option<StateVecEdits<T>>,
}

impl<T> IdxAdjustment<T> {
    /// Returns an [`IdxAdjustment`] that advances to the next index without
    /// removing any states.
    #[inline]
    #[must_use]
    pub const fn next() -> Self {
        Self {
            advance: true,
            edits:   None,
        }
    }
}

/// Flags indicating which of the states in [`StateWithFlanking`] should be
/// removed, after the modifications made by the closure in [`rewrite`].
///
/// Note that no `remove_left2` field is present. `remove_right2` is provided in
/// case `remove_right1` is removed (e.g., due to shrinking to 0) and then
/// `right2` needs to be merged into `current` (and hence be removed). When the
/// same happens but on the left, `current` is merged into `left2`, and then
/// `current` is removed. Hence, there is not typically a need for removing
/// `left2`.
#[derive(Clone, Eq, PartialEq, Debug)]
pub struct StateVecEdits<T> {
    /// The state left of the current one (`idx-1`) should be removed.
    pub remove_left1:   bool,
    /// The current state (`idx`) should be removed.
    pub remove_current: bool,
    /// The state right of the current one (`idx+1`) should be removed.
    pub remove_right1:  bool,
    /// The state two to the right of the current one (`idx+2`) should be
    /// removed.
    pub remove_right2:  bool,
    /// After removal, a state to insert left of the current one (before the
    /// value at `idx`).
    pub insert_left:    Option<T>,
    /// After removal, a state to insert right of the current one (after the
    /// value at `idx`).
    pub insert_right:   Option<T>,
}

impl<T> Default for StateVecEdits<T> {
    fn default() -> Self {
        Self {
            remove_left1:   false,
            remove_current: false,
            remove_right1:  false,
            remove_right2:  false,
            insert_left:    None,
            insert_right:   None,
        }
    }
}

impl<T> StateVecEdits<T> {
    /// The return value is the amount that must be added to `idx` to correct
    /// for the removed/inserted states.
    ///
    /// The shift will alter the index so that it corresponds to the index
    /// before the next state to be processed. The next state to be processed is
    /// the leftmost non-removed state starting at `right1`.
    ///
    /// For example, if no insertions occur and `current` is not removed, the
    /// shift will cause the index to remain at `current`. If `current` is
    /// removed and no insertions occur, then index will correspond to `left1`,
    /// or `left2` is `left1` was removed. If `current` is not removed and a
    /// state is inserted on the left, then the index remains at `current`. If
    /// the insertion is instead on the right, then the index is now at the
    /// newly inserted right state.
    ///
    /// In this way, assuming the index gets advanced by the caller, neither
    /// `current` nor any inserted states get reprocessed.
    pub fn apply(self, states: &mut Vec<T>, idx: usize) -> isize {
        let StateVecEdits {
            remove_left1,
            remove_current,
            remove_right1,
            remove_right2,
            insert_left,
            insert_right,
        } = self;

        let mut idx_offset = 0;

        // Alter states from right to left
        if remove_right2 {
            states.remove(idx + 2);
        }
        if remove_right1 {
            states.remove(idx + 1);
        }
        if let Some(insert_right) = insert_right {
            states.insert(idx + 1, insert_right);
            idx_offset += 1;
        }
        if remove_current {
            states.remove(idx);
            idx_offset -= 1;
        }
        if let Some(insert_left) = insert_left {
            states.insert(idx, insert_left);
            idx_offset += 1;
        }
        if remove_left1 {
            states.remove(idx - 1);
            idx_offset -= 1;
        }

        idx_offset
    }
}
