# Changelog

All notable changes to `lambert_izzo` are documented here. The format
follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and
the project adheres to [Semantic Versioning](https://semver.org/) once
it reaches `1.0`.

## [Unreleased]

Targeting a `1.0.0` release once the API and feature surface have settled
through external review.

## [0.4.0] — 2026-04-26

### Added

- **`no_std` support.** Core kernel compiles without the standard
  library; transcendental math routes through `num_traits::Float` (with
  `libm`). Tested via `cargo build --target wasm32-unknown-unknown -p
  lambert_izzo --no-default-features --lib`.
- **`serde` feature.** Optional `Serialize`/`Deserialize` derives on
  every public type, including `LambertError` (now serializes to a
  discriminated union via the externally-tagged default).
- **`test-utils` feature.** Promotes the universal-variable Kepler
  propagator from `#[cfg(test)] pub(crate)` to a public function under
  `lambert_izzo::test_utils::kepler_propagate` so downstream
  integration tests can round-trip-validate Lambert solutions without
  re-implementing it.
- **`rayon` feature.** Adds `lambert_par_iter` for parallel batch
  evaluation. Incompatible with `no_std` (Rayon pulls `std`).
- **Batch streaming API.** New `LambertInput` struct and
  `lambert_iter(&[LambertInput])` for porkchop-plot-style workloads.
  Allocation-free.
- **Criterion benchmarks** under `crates/lambert_izzo/benches/`:
  `single_rev`, `multi_rev`, `batch` (with `rayon` toggle).
- **MSRV declared.** `rust-version = "1.85"` (the first edition-2024
  stable release). Lowered from the conservative `1.88` pin.

### Changed

- **`LambertError::NonFiniteInput::parameter`** field type changes from
  `&'static str` to a typed `NonFiniteParameter` enum. Required for
  serde round-tripping under `no_std`; also enables exhaustive matching
  in callers.

### WASM adapter (`lambert_izzo_wasm` v0.3.0)

- **Structured errors.** `solveLambert(...)` now rejects with a typed
  discriminated union (`LambertErrorOutput`, `Tsify`-derived) instead
  of a plain string. JS callers can `switch` on `kind` without parsing
  error messages.

## [0.3.0] — 2026-04-26

### Changed (breaking)

- **Public surface uses `[f64; 3]`** for all position/velocity vectors.
  The hard `nalgebra` dependency is dropped from the runtime;
  `nalgebra::Vector3<f64>` and `glam::DVec3` already convert to/from
  `[f64; 3]` natively.
- **`Vec<LambertSolution>` → `LambertSolutions { single, multi:
  ArrayVec<MultiRevPair, MAX_MULTI_REV_PAIRS> }`.** The implicit
  chunk-by-2 multi-rev convention is gone; pairing is encoded in the
  `MultiRevPair { n_revs, long_period, short_period }` type. Zero heap
  allocation on the solver hot path.
- **`SolverDiagnostics` removed from `LambertSolution`.** A new
  `solve_with_diagnostics(...)` returns `(LambertSolutions,
  LambertDiagnostics)` for callers that need iteration counts and the
  Lancaster–Blanchard `x`.

### Added

- `lambert_both_ways(...)` — single call returning both transfer
  directions, for the porkchop-plot pattern.

### WASM adapter (`lambert_izzo_wasm` v0.2.0)

- Array-conversion shim deleted; the public surface was already arrays.

## [0.2.0] — 2026-04-26

- Workspace restructure: `crates/lambert_izzo` (core) and
  `crates/lambert_izzo_wasm` (`wasm-bindgen` adapter).
- Documentation polish on `pub fn lambert`; tightened `Cargo.toml`
  metadata for `crates.io` publication.

## [0.1.0] — 2026-04-26

- Initial implementation of Izzo's revisited Lambert solver
  (single + multi-revolution, short + long way).
- Three-regime time-of-flight dispatch (Battin / Lancaster–Blanchard /
  Lagrange) with documented thresholds in `constants.rs`.
- Householder iteration with derivative-matched starters (Izzo Eq. 30,
  31) and Halley `T_min` search for multi-rev infeasibility.
- Round-trip Kepler validation via universal-variable propagator;
  `examples/stress.rs` reproduces Izzo §5 statistical sweeps.
- Strict lint baseline: `clippy::pedantic` + bans on `unwrap`/`expect`/
  `panic`/`unreachable` in lib code.

[Unreleased]: https://github.com/sakobu/izzos_lambert/compare/v0.4.0...HEAD
[0.4.0]: https://github.com/sakobu/izzos_lambert/compare/v0.3.0...v0.4.0
[0.3.0]: https://github.com/sakobu/izzos_lambert/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/sakobu/izzos_lambert/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/sakobu/izzos_lambert/releases/tag/v0.1.0
