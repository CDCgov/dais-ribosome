use crate::{
    QueryRecord,
    annotation::rewrites::{
        fix_frames::{ShiftDir, pick_del_shift_with_stats, pick_ins_shift_with_stats},
        get_states::StateVecEdits,
    },
    config::ProductSpec,
    data::{exons::Exons, weights::CodonPositionWeights},
    ranges::{CdsDeletionRange, CdsInsertionRange, CdsMatchRange, InsertionIdx},
};

/// A purely empty [`Exons`] struct for use in the below unit-tests. Note that
/// this invalidates the assumptions on [`Exons`], but is acceptable for these
/// tests since the `exons` field is not accessed.
const EMPTY_EXONS: Exons = Exons {
    required_start:     None,
    coords:             Vec::new(),
    overlapped_regions: Vec::new(),
    noncoding_regions:  Vec::new(),
};

#[test]
fn shift_ins_left1() {
    let mut codon_weights = CodonPositionWeights::new();
    codon_weights.insert(1, *b"ACG", 1).unwrap();
    codon_weights.insert(1, *b"TAC", 10).unwrap();
    let product_spec = ProductSpec {
        name:          String::from("shift_ins_left1"),
        exons:         EMPTY_EXONS,
        codon_weights: Some(codon_weights),
    };

    let left_match = CdsMatchRange {
        query_range: 0..1,
        cds_range:   0..1,
    };

    let ins = CdsInsertionRange {
        cds_index:   InsertionIdx::from_left_idx(0),
        query_range: 1..4,
    };

    let right_match = CdsMatchRange {
        query_range: 4..6,
        cds_range:   1..3,
    };

    let query = QueryRecord::new(String::from("test_seq"), b"ACGTAC".to_vec(), String::from("test_ctype")).unwrap();

    let dir = pick_ins_shift_with_stats(&left_match, &ins, &right_match, &query, &product_spec);

    assert_eq!(dir, Some(ShiftDir::Left))
}

#[test]
fn shift_ins_left1_default() {
    let mut codon_weights = CodonPositionWeights::new();
    codon_weights.insert(1, *b"ACG", 1).unwrap();
    codon_weights.insert(1, *b"TAC", 1).unwrap();
    let product_spec = ProductSpec {
        name:          String::from("shift_ins_left1_default"),
        exons:         EMPTY_EXONS,
        codon_weights: Some(codon_weights),
    };

    let left_match = CdsMatchRange {
        query_range: 0..1,
        cds_range:   0..1,
    };

    let ins = CdsInsertionRange {
        cds_index:   InsertionIdx::from_left_idx(0),
        query_range: 1..4,
    };

    let right_match = CdsMatchRange {
        query_range: 4..6,
        cds_range:   1..3,
    };

    let query = QueryRecord::new(String::from("test_seq"), b"ACGTAC".to_vec(), String::from("test_ctype")).unwrap();

    let dir = pick_ins_shift_with_stats(&left_match, &ins, &right_match, &query, &product_spec);

    assert_eq!(dir, Some(ShiftDir::Left))
}

#[test]
fn shift_ins_left2() {
    let mut codon_weights = CodonPositionWeights::new();
    codon_weights.insert(1, *b"ACG", 1).unwrap();
    codon_weights.insert(1, *b"TAC", 10).unwrap();
    let product_spec = ProductSpec {
        name:          String::from("shift_ins_left2"),
        exons:         EMPTY_EXONS,
        codon_weights: Some(codon_weights),
    };

    let left_match = CdsMatchRange {
        query_range: 0..2,
        cds_range:   0..2,
    };

    let ins = CdsInsertionRange {
        cds_index:   InsertionIdx::from_left_idx(1),
        query_range: 2..5,
    };

    let right_match = CdsMatchRange {
        query_range: 5..6,
        cds_range:   2..3,
    };

    let query = QueryRecord::new(String::from("test_seq"), b"ACGTAC".to_vec(), String::from("test_ctype")).unwrap();

    let dir = pick_ins_shift_with_stats(&left_match, &ins, &right_match, &query, &product_spec);

    assert_eq!(dir, Some(ShiftDir::Left))
}

#[test]
fn shift_ins_right1() {
    let mut codon_weights = CodonPositionWeights::new();
    codon_weights.insert(1, *b"ACG", 10).unwrap();
    codon_weights.insert(1, *b"TAC", 1).unwrap();
    let product_spec = ProductSpec {
        name:          String::from("shift_ins_right1"),
        exons:         EMPTY_EXONS,
        codon_weights: Some(codon_weights),
    };

    let left_match = CdsMatchRange {
        query_range: 0..2,
        cds_range:   0..2,
    };

    let ins = CdsInsertionRange {
        cds_index:   InsertionIdx::from_left_idx(1),
        query_range: 2..5,
    };

    let right_match = CdsMatchRange {
        query_range: 5..6,
        cds_range:   2..3,
    };

    let query = QueryRecord::new(String::from("test_seq"), b"ACGTAC".to_vec(), String::from("test_ctype")).unwrap();

    let dir = pick_ins_shift_with_stats(&left_match, &ins, &right_match, &query, &product_spec);

    assert_eq!(dir, Some(ShiftDir::Right))
}

#[test]
fn shift_ins_right1_default() {
    let mut codon_weights = CodonPositionWeights::new();
    codon_weights.insert(1, *b"ACG", 1).unwrap();
    codon_weights.insert(1, *b"TAC", 1).unwrap();
    let product_spec = ProductSpec {
        name:          String::from("shift_ins_right1_default"),
        exons:         EMPTY_EXONS,
        codon_weights: Some(codon_weights),
    };

    let left_match = CdsMatchRange {
        query_range: 0..2,
        cds_range:   0..2,
    };

    let ins = CdsInsertionRange {
        cds_index:   InsertionIdx::from_left_idx(1),
        query_range: 2..5,
    };

    let right_match = CdsMatchRange {
        query_range: 5..6,
        cds_range:   2..3,
    };

    let query = QueryRecord::new(String::from("test_seq"), b"ACGTAC".to_vec(), String::from("test_ctype")).unwrap();

    let dir = pick_ins_shift_with_stats(&left_match, &ins, &right_match, &query, &product_spec);

    assert_eq!(dir, Some(ShiftDir::Right))
}

#[test]
fn shift_ins_right2() {
    let mut codon_weights = CodonPositionWeights::new();
    codon_weights.insert(1, *b"ACG", 10).unwrap();
    codon_weights.insert(1, *b"TAC", 1).unwrap();
    let product_spec = ProductSpec {
        name:          String::from("shift_ins_right2"),
        exons:         EMPTY_EXONS,
        codon_weights: Some(codon_weights),
    };

    let left_match = CdsMatchRange {
        query_range: 0..1,
        cds_range:   0..1,
    };

    let ins = CdsInsertionRange {
        cds_index:   InsertionIdx::from_left_idx(0),
        query_range: 1..4,
    };

    let right_match = CdsMatchRange {
        query_range: 4..6,
        cds_range:   1..3,
    };

    let query = QueryRecord::new(String::from("test_seq"), b"ACGTAC".to_vec(), String::from("test_ctype")).unwrap();

    let dir = pick_ins_shift_with_stats(&left_match, &ins, &right_match, &query, &product_spec);

    assert_eq!(dir, Some(ShiftDir::Right))
}

#[test]
fn shift_del_left1() {
    let mut codon_weights = CodonPositionWeights::new();
    codon_weights.insert(1, *b"ACG", 1).unwrap();
    codon_weights.insert(2, *b"ACG", 10).unwrap();
    let product_spec = ProductSpec {
        name:          String::from("shift_del_left1"),
        exons:         EMPTY_EXONS,
        codon_weights: Some(codon_weights),
    };

    let left_match = CdsMatchRange {
        query_range: 0..1,
        cds_range:   0..1,
    };

    let del = CdsDeletionRange { cds_range: 1..4 };

    let right_match = CdsMatchRange {
        query_range: 1..3,
        cds_range:   4..6,
    };

    let query = QueryRecord::new(String::from("test_seq"), b"ACG".to_vec(), String::from("test_ctype")).unwrap();

    let dir = pick_del_shift_with_stats(&left_match, &del, &right_match, &query, &product_spec);

    assert_eq!(dir, Some(ShiftDir::Left))
}

#[test]
fn shift_del_left1_default() {
    let mut codon_weights = CodonPositionWeights::new();
    codon_weights.insert(1, *b"ACG", 1).unwrap();
    codon_weights.insert(2, *b"ACG", 1).unwrap();
    let product_spec = ProductSpec {
        name:          String::from("shift_del_left1_default"),
        exons:         EMPTY_EXONS,
        codon_weights: Some(codon_weights),
    };

    let left_match = CdsMatchRange {
        query_range: 0..1,
        cds_range:   0..1,
    };

    let del = CdsDeletionRange { cds_range: 1..4 };

    let right_match = CdsMatchRange {
        query_range: 1..3,
        cds_range:   4..6,
    };

    let query = QueryRecord::new(String::from("test_seq"), b"ACG".to_vec(), String::from("test_ctype")).unwrap();

    let dir = pick_del_shift_with_stats(&left_match, &del, &right_match, &query, &product_spec);

    assert_eq!(dir, Some(ShiftDir::Left))
}

#[test]
fn shift_del_left2() {
    let mut codon_weights = CodonPositionWeights::new();
    codon_weights.insert(1, *b"ACG", 1).unwrap();
    codon_weights.insert(2, *b"ACG", 10).unwrap();
    let product_spec = ProductSpec {
        name:          String::from("shift_del_left2"),
        exons:         EMPTY_EXONS,
        codon_weights: Some(codon_weights),
    };

    let left_match = CdsMatchRange {
        query_range: 0..2,
        cds_range:   0..2,
    };

    let del = CdsDeletionRange { cds_range: 2..5 };

    let right_match = CdsMatchRange {
        query_range: 2..3,
        cds_range:   5..6,
    };

    let query = QueryRecord::new(String::from("test_seq"), b"ACG".to_vec(), String::from("test_ctype")).unwrap();

    let dir = pick_del_shift_with_stats(&left_match, &del, &right_match, &query, &product_spec);

    assert_eq!(dir, Some(ShiftDir::Left))
}

#[test]
fn shift_del_right1() {
    let mut codon_weights = CodonPositionWeights::new();
    codon_weights.insert(1, *b"ACG", 10).unwrap();
    codon_weights.insert(2, *b"ACG", 1).unwrap();
    let product_spec = ProductSpec {
        name:          String::from("shift_del_right1"),
        exons:         EMPTY_EXONS,
        codon_weights: Some(codon_weights),
    };

    let left_match = CdsMatchRange {
        query_range: 0..2,
        cds_range:   0..2,
    };

    let del = CdsDeletionRange { cds_range: 2..5 };

    let right_match = CdsMatchRange {
        query_range: 2..3,
        cds_range:   5..6,
    };

    let query = QueryRecord::new(String::from("test_seq"), b"ACG".to_vec(), String::from("test_ctype")).unwrap();

    let dir = pick_del_shift_with_stats(&left_match, &del, &right_match, &query, &product_spec);

    assert_eq!(dir, Some(ShiftDir::Right))
}

#[test]
fn shift_del_right1_default() {
    let mut codon_weights = CodonPositionWeights::new();
    codon_weights.insert(1, *b"ACG", 1).unwrap();
    codon_weights.insert(2, *b"ACG", 1).unwrap();
    let product_spec = ProductSpec {
        name:          String::from("shift_del_right1_default"),
        exons:         EMPTY_EXONS,
        codon_weights: Some(codon_weights),
    };

    let left_match = CdsMatchRange {
        query_range: 0..2,
        cds_range:   0..2,
    };

    let del = CdsDeletionRange { cds_range: 2..5 };

    let right_match = CdsMatchRange {
        query_range: 2..3,
        cds_range:   5..6,
    };

    let query = QueryRecord::new(String::from("test_seq"), b"ACG".to_vec(), String::from("test_ctype")).unwrap();

    let dir = pick_del_shift_with_stats(&left_match, &del, &right_match, &query, &product_spec);

    assert_eq!(dir, Some(ShiftDir::Right))
}

#[test]
fn shift_del_right2() {
    let mut codon_weights = CodonPositionWeights::new();
    codon_weights.insert(1, *b"ACG", 10).unwrap();
    codon_weights.insert(2, *b"ACG", 1).unwrap();
    let product_spec = ProductSpec {
        name:          String::from("shift_del_right2"),
        exons:         EMPTY_EXONS,
        codon_weights: Some(codon_weights),
    };

    let left_match = CdsMatchRange {
        query_range: 0..1,
        cds_range:   0..1,
    };

    let del = CdsDeletionRange { cds_range: 1..4 };

    let right_match = CdsMatchRange {
        query_range: 1..3,
        cds_range:   4..6,
    };

    let query = QueryRecord::new(String::from("test_seq"), b"ACG".to_vec(), String::from("test_ctype")).unwrap();

    let dir = pick_del_shift_with_stats(&left_match, &del, &right_match, &query, &product_spec);

    assert_eq!(dir, Some(ShiftDir::Right))
}

#[test]
fn state_vec_edit_no_change() {
    let edits = StateVecEdits::default();
    for idx in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, usize::MAX] {
        let mut states = vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9];
        let shift = edits.clone().apply(&mut states, idx);
        assert_eq!(shift, 0);
        assert_eq!(states, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]);
    }
}

#[test]
fn state_vec_edit() {
    /// Abstracts the logic for this unit test. The states being edited are 0-9.
    /// Inserts on left should be 10, and inserts on right should be 11. Edits
    /// are applied at index 3.
    fn return_applied(edits: StateVecEdits<u8>) -> Vec<u8> {
        let mut states = vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9];
        let idx = 3;

        let shift = edits.apply(&mut states, idx);
        let new_idx = states.iter().position(|val| (idx as u8 + 1..=9).contains(val)).unwrap() - 1;
        assert_eq!(idx.wrapping_add_signed(shift), new_idx);

        states
    }

    // Remove L1
    let edits = StateVecEdits {
        remove_left1: true,
        ..Default::default()
    };
    let expected = [0, 1, 3, 4, 5, 6, 7, 8, 9];
    assert_eq!(return_applied(edits), expected);

    // Remove C
    let edits = StateVecEdits {
        remove_current: true,
        ..Default::default()
    };
    let expected = [0, 1, 2, 4, 5, 6, 7, 8, 9];
    assert_eq!(return_applied(edits), expected);

    // Remove R1
    let edits = StateVecEdits {
        remove_right1: true,
        ..Default::default()
    };
    let expected = [0, 1, 2, 3, 5, 6, 7, 8, 9];
    assert_eq!(return_applied(edits), expected);

    // Remove R2
    let edits = StateVecEdits {
        remove_right2: true,
        ..Default::default()
    };
    let expected = [0, 1, 2, 3, 4, 6, 7, 8, 9];
    assert_eq!(return_applied(edits), expected);

    // Insert left
    let edits = StateVecEdits {
        insert_left: Some(10),
        ..Default::default()
    };
    let expected = [0, 1, 2, 10, 3, 4, 5, 6, 7, 8, 9];
    assert_eq!(return_applied(edits), expected);

    // Insert right
    let edits = StateVecEdits {
        insert_right: Some(11),
        ..Default::default()
    };
    let expected = [0, 1, 2, 3, 11, 4, 5, 6, 7, 8, 9];
    assert_eq!(return_applied(edits), expected);

    // Insert both
    let edits = StateVecEdits {
        insert_left: Some(10),
        insert_right: Some(11),
        ..Default::default()
    };
    let expected = [0, 1, 2, 10, 3, 11, 4, 5, 6, 7, 8, 9];
    assert_eq!(return_applied(edits), expected);

    // Remove C, insert left
    let edits = StateVecEdits {
        remove_current: true,
        insert_left: Some(10),
        ..Default::default()
    };
    let expected = [0, 1, 2, 10, 4, 5, 6, 7, 8, 9];
    assert_eq!(return_applied(edits), expected);

    // Remove C, insert right
    let edits = StateVecEdits {
        remove_current: true,
        insert_right: Some(11),
        ..Default::default()
    };
    let expected = [0, 1, 2, 11, 4, 5, 6, 7, 8, 9];
    assert_eq!(return_applied(edits), expected);

    // Remove C, insert both
    let edits = StateVecEdits {
        remove_current: true,
        insert_left: Some(10),
        insert_right: Some(11),
        ..Default::default()
    };
    let expected = [0, 1, 2, 10, 11, 4, 5, 6, 7, 8, 9];
    assert_eq!(return_applied(edits), expected);

    // Remove L1, insert left
    let edits = StateVecEdits {
        remove_left1: true,
        insert_left: Some(10),
        ..Default::default()
    };
    let expected = [0, 1, 10, 3, 4, 5, 6, 7, 8, 9];
    assert_eq!(return_applied(edits), expected);

    // Remove L1, insert right
    let edits = StateVecEdits {
        remove_left1: true,
        insert_right: Some(11),
        ..Default::default()
    };
    let expected = [0, 1, 3, 11, 4, 5, 6, 7, 8, 9];
    assert_eq!(return_applied(edits), expected);

    // Remove L1, insert both
    let edits = StateVecEdits {
        remove_left1: true,
        insert_left: Some(10),
        insert_right: Some(11),
        ..Default::default()
    };
    let expected = [0, 1, 10, 3, 11, 4, 5, 6, 7, 8, 9];
    assert_eq!(return_applied(edits), expected);

    // Remove R1, insert left
    let edits = StateVecEdits {
        remove_right1: true,
        insert_left: Some(10),
        ..Default::default()
    };
    let expected = [0, 1, 2, 10, 3, 5, 6, 7, 8, 9];
    assert_eq!(return_applied(edits), expected);

    // Remove R1, insert right
    let edits = StateVecEdits {
        remove_right1: true,
        insert_right: Some(11),
        ..Default::default()
    };
    let expected = [0, 1, 2, 3, 11, 5, 6, 7, 8, 9];
    assert_eq!(return_applied(edits), expected);

    // Remove R1, insert both
    let edits = StateVecEdits {
        remove_right1: true,
        insert_left: Some(10),
        insert_right: Some(11),
        ..Default::default()
    };
    let expected = [0, 1, 2, 10, 3, 11, 5, 6, 7, 8, 9];
    assert_eq!(return_applied(edits), expected);

    // Remove R2, insert left
    let edits = StateVecEdits {
        remove_right2: true,
        insert_left: Some(10),
        ..Default::default()
    };
    let expected = [0, 1, 2, 10, 3, 4, 6, 7, 8, 9];
    assert_eq!(return_applied(edits), expected);

    // Remove R2, insert right
    let edits = StateVecEdits {
        remove_right2: true,
        insert_right: Some(11),
        ..Default::default()
    };
    let expected = [0, 1, 2, 3, 11, 4, 6, 7, 8, 9];
    assert_eq!(return_applied(edits), expected);

    // Remove R2, insert both
    let edits = StateVecEdits {
        remove_right2: true,
        insert_left: Some(10),
        insert_right: Some(11),
        ..Default::default()
    };
    let expected = [0, 1, 2, 10, 3, 11, 4, 6, 7, 8, 9];
    assert_eq!(return_applied(edits), expected);

    // Remove L1 and R1, insert left
    let edits = StateVecEdits {
        remove_left1: true,
        remove_right1: true,
        insert_left: Some(10),
        ..Default::default()
    };
    let expected = [0, 1, 10, 3, 5, 6, 7, 8, 9];
    assert_eq!(return_applied(edits), expected);

    // Remove L1 and R1, insert right
    let edits = StateVecEdits {
        remove_left1: true,
        remove_right1: true,
        insert_right: Some(11),
        ..Default::default()
    };
    let expected = [0, 1, 3, 11, 5, 6, 7, 8, 9];
    assert_eq!(return_applied(edits), expected);

    // Remove L1 and R1, insert both
    let edits = StateVecEdits {
        remove_left1: true,
        remove_right1: true,
        insert_left: Some(10),
        insert_right: Some(11),
        ..Default::default()
    };
    let expected = [0, 1, 10, 3, 11, 5, 6, 7, 8, 9];
    assert_eq!(return_applied(edits), expected);

    // Remove L1 and R2, insert left
    let edits = StateVecEdits {
        remove_left1: true,
        remove_right2: true,
        insert_left: Some(10),
        ..Default::default()
    };
    let expected = [0, 1, 10, 3, 4, 6, 7, 8, 9];
    assert_eq!(return_applied(edits), expected);

    // Remove L1 and R2, insert right
    let edits = StateVecEdits {
        remove_left1: true,
        remove_right2: true,
        insert_right: Some(11),
        ..Default::default()
    };
    let expected = [0, 1, 3, 11, 4, 6, 7, 8, 9];
    assert_eq!(return_applied(edits), expected);

    // Remove L1 and R2, insert both
    let edits = StateVecEdits {
        remove_left1: true,
        remove_right2: true,
        insert_left: Some(10),
        insert_right: Some(11),
        ..Default::default()
    };
    let expected = [0, 1, 10, 3, 11, 4, 6, 7, 8, 9];
    assert_eq!(return_applied(edits), expected);

    // Remove C and L1, insert left
    let edits = StateVecEdits {
        remove_left1: true,
        remove_current: true,
        insert_left: Some(10),
        ..Default::default()
    };
    let expected = [0, 1, 10, 4, 5, 6, 7, 8, 9];
    assert_eq!(return_applied(edits), expected);

    // Remove C and L1, insert right
    let edits = StateVecEdits {
        remove_left1: true,
        remove_current: true,
        insert_right: Some(11),
        ..Default::default()
    };
    let expected = [0, 1, 11, 4, 5, 6, 7, 8, 9];
    assert_eq!(return_applied(edits), expected);

    // Remove C and L1, insert both
    let edits = StateVecEdits {
        remove_left1: true,
        remove_current: true,
        insert_left: Some(10),
        insert_right: Some(11),
        ..Default::default()
    };
    let expected = [0, 1, 10, 11, 4, 5, 6, 7, 8, 9];
    assert_eq!(return_applied(edits), expected);

    // Remove C and R1, insert left
    let edits = StateVecEdits {
        remove_current: true,
        remove_right1: true,
        insert_left: Some(10),
        ..Default::default()
    };
    let expected = [0, 1, 2, 10, 5, 6, 7, 8, 9];
    assert_eq!(return_applied(edits), expected);

    // Remove C and R1, insert right
    let edits = StateVecEdits {
        remove_current: true,
        remove_right1: true,
        insert_right: Some(11),
        ..Default::default()
    };
    let expected = [0, 1, 2, 11, 5, 6, 7, 8, 9];
    assert_eq!(return_applied(edits), expected);

    // Remove C and R1, insert both
    let edits = StateVecEdits {
        remove_current: true,
        remove_right1: true,
        insert_left: Some(10),
        insert_right: Some(11),
        ..Default::default()
    };
    let expected = [0, 1, 2, 10, 11, 5, 6, 7, 8, 9];
    assert_eq!(return_applied(edits), expected);
}
