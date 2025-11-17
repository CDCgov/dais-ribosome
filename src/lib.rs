#![feature(try_trait_v2)]
#![allow(dead_code)]

use rand::{Rng, distr::Alphanumeric};
use zoe::prelude::*;

// A well-formed DAIS-ribosome input.
struct InputRecord {
    id:    String,
    taxon: String,
    seq:   Nucleotides,
}

pub fn generate_token(length: usize) -> String {
    rand::rng().sample_iter(&Alphanumeric).take(length).map(char::from).collect()
}

pub mod conf;
