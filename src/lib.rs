#![feature(try_trait_v2, int_format_into)]

// =============================================================================
// Module hierarchy
// =============================================================================

/// Configuration loading and path resolution.
pub mod config;

/// Data loading and structures for module resources.
pub mod data;

/// Protein annotation engine.
pub mod annotation;
