//! Outputs from DAIS-ribosome annotation, including range-based structs and
//! materialized/computed types.

mod computed_genome;
mod computed_product;
mod output;

pub use computed_genome::*;
pub use computed_product::*;
pub use output::*;
