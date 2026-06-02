use crate::{
    QueryRecord,
    annotation::fix_frames::{ShiftDir, pick_del_shift_with_stats, pick_ins_shift_with_stats},
    config::ProductSpec,
    data::{exons::Exons, weights::CodonPositionWeights},
    ranges::{CdsDeletionRange, CdsInsertionRange, CdsMatchRange, InsertionIdx},
};

/// A purely empty [`Exons`] struct for use in the below unit-tests. Note
/// that this invalidates the assumptions on [`Exons`], but is acceptable
/// for these tests since the `exons` field is not accessed.
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
