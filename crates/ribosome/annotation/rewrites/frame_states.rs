/// The state at a given index, along with flanking states if they are in
/// bounds.
///
/// This is returned by [`get_frame_states`], and used by [`fix_frame`].
///
/// ## Parameters
///
/// - `'a`: The lifetime of the mutable reference
/// - `T`: The type of the range, such as [`CdsStateRange`] or [`StateRange`]
pub struct FrameStates<'a, T> {
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
pub fn get_frame_states<T>(idx: usize, product_ranges: &mut [T]) -> Option<FrameStates<'_, T>> {
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

    Some(FrameStates {
        left2,
        left1,
        current,
        right1,
        right2,
    })
}
