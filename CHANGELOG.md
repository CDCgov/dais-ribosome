# DAIS-ribosome Changelog

All notable changes to this project will be documented in this file. The format
is roughly based on [Keep a Changelog], and this project adheres to
[Semantic Versioning].

## [2.0.0] - 2026-07-08

### Changed

- DAIS-ribosome has been rewritten into self-contained Rust CLI app. See [`MIGRATION.md`] for changes between v1 and v2.

## [1.7.1] - 2026-05-14

### Fixed

- Fixes a bug where some insertions after the first may not be spliced into the unaligned sequence properly
- The translation of insertions now supports `U`
- An out-of-frame deletion that occurs within 3 bases of another is properly shifted
- Partially fixes a bug where trailing insertions appear at the end of a product. **It is a known bug that in some cases, an insertion still may appear after a stop codon in the aligned CDS sequence for a product. This will be resolved in DAIS-ribosome v2**
- Trailing insertions at the end of a product are now properly excluded from the
  unaligned sequence outputs (fix in `convert`)

## [1.7.0] - 2026-04-23

### Added

- DAIS-ribosome now supports aarch64/arm64 on both Mac and Linux, albeit SSWSORT requires Rosetta 2 if it is called.

### Changed

- The ranges for the `HA1` product for `A_HA_H1` and `A_HA_H3` are extended by 3 nucleotides and 18 nucleotides respectively
- In advance of ribosome/sswsort v2, we pin ribosome v1 to sswsort v1.
- Codons containing a mix of `.` and `-` now translate to `.` instead of `-` (fix in `convert`)

### Fixed

- The check for the required start codon in PB1-F2 now allows for AUG in addition to ATG
- Fixes a bug causing some query coordinate ranges to be incorrectly shifted right for `M2`, `NEP`, and `PA-X`
- Fixes Dockerfile build process to allow for custom certificates
- Fixes the majority of bugs where indels in non-coding regions are not reflected in query coordinates. **It is a known bug that in rare cases, the query coordinates may still be incorrect. This will be resolved in DAIS-ribosome v2**
- Trailing indels that are then followed by padding no longer are counted as
  causing a frame shift (fix in `convert`)
- Stop extension now properly occurs even in the presence of prior insertions (fix in `convert`)
- Insertions near the end of a product are now properly output, even if their
  frames could not be corrected (fix in `editMSA`)
- The default codon/position weights are now properly accessed when shifting
  insertions, which may improve the direction chosen when position-specific
  information is unavailable or has a tie (fix in `editMSA`)
- Insertions that shift to the beginning of a product are now properly removed (fix in `editMSA`)
- Deletions that are near each other are now both properly shifted, instead of
  just the first (fix in `editMSA`)

## [1.6.2] - 2026-03-19

### Fixed

- The check for the required start codon in PB1-F2 properly checks the beginning of the coding sequence instead of the query sequence
- Translation from nucleotides to amino acids for matched regions properly supports `U` (fix in `convert`)

## [1.6.1] - 2025-04-18

### Fixed

- The codon-weight matrix had the wrong number of sites for KF640637|L, this has been corrected.

## [1.6.0] - 2025-04-18

### Added

- Adds the A_NA_N4 reference, 'ONTARIO6118' to the INFLUENZA module.
- Adds the A_NA_N5 reference, 'RU1526' and 'ALASKA4733' (different lineage coordinate spaces) to the INFLUENZA module.
- Adds the A_NA_N6 reference, 'SICHUAN26221' to the INFLUENZA module.

## [1.5.6] - 2025-03-28

### Changed

- ASTRAKHAN3212|A_NA_N8|OM403994 references now use lineage representatives LC339685 and MF973227, which enhances alignment quality.

## [1.5.5] - 2024-10-30

### Changed

- Streamlines the dockerfile and changes to debian/bookworm.

## [1.5.4] - 2024-08-22

### Changed

- RSV ctypes are renamed, e.g., RSVA is RSV_AD; RSV reference IDs now match the NCBI accession.

### Fixed

- Updates the SC2 codon-weight matrix to induce corrected alignments for a deletion in S circa position 24.

## [1.5.3] - 2024-07-25

### Fixed

- SC2 alignments were sometimes missing biologically relevant indels. Smith-Waterman weights were relaxed to accommodate.

## [1.5.2] - 2024-05-10

### Added

- Adds docker support and adds usage.
- Adds configurable work directories using `IFX_WORK_DIR`.

### Fixed

- Fixes SC2 alignments for S protein insertions (MPLF) at position 16.

## [1.5.1] - 2024-04-30

### Fixed

- RSV annotations now correctly validated.

## [1.5.0] - 2024-04-25

### Added

- Adds the A_NA_N8 reference, 'ASTRAKHAN3212' to the INFLUENZA module.

## [1.4.0] - 2024-04-03

### Added

- Adds built-in RSV module support. Many Thanks to C. Paden and P. Mandal!

## [1.3.3] - 2023-06-26

### Changed

- Right pads the CDS after premature stop codons.

### Fixed

- Corrects coordinate map / range edge case involving insertions and incomplete CDS.
- "has insertion" boolean now reports `false` when the only insertions are outside the terminated CDS.
- The CDS insertion table no longer provides insertions outside of the terminated CDS.

## [1.3.2] - 2023-06-23

### Changed

- Relaxes input file format detection and makes it invariant for all modules.

## [1.3.1] - 2023-05-24

### Changed

- Restores right padding of CDS without adding back in the variable stop codon bugs from v1.2.
- Premature stop codons will still be truncated.

## [1.3.0] - 2023-05

### Fixed

- Parallel SSW errors bubble up properly (including for fail-over from GRID jobs).
- Properly terminates CDS similar to AA sequences for premature stop codons.
- Skip empty exons and sequences during insertion processing.

## [1.2.0] - 2022-08

### Added

- Adds ORF7b and ORF9b to the BETACORONAVIRUS.spec (Thanks to C. Paden).

## [1.1.0] - 2022-07

### Changed

- Cleaned up scripts and added individual licensing.

## [1.0.0] - 2021-04

### Added

- Initial tagged release supporting SARS-CoV-2 and Influenza.

<!-- Versions -->

[2.0.0]: https://github.com/CDCgov/dais-ribosome/compare/v1.7.1...v2.0.0
[keep a changelog]: https://keepachangelog.com/en/1.0.0/
[semantic versioning]: https://semver.org/spec/v2.0.0.html
