# Changelog

All notable changes to `lambert_izzo` are documented here. The format
follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and
the project adheres to [Semantic Versioning](https://semver.org/) once
it reaches `1.0`.

## [Unreleased]

Targeting a `1.0.0` release once the API and feature surface have settled
through external review.

### Removed (breaking)

- **`test-utils` feature** removed from `lambert_izzo`. The Kepler
  propagator it exposed has moved to a new workspace-internal crate,
  `lambert_izzo_test_support` (path-only, `publish = false`). Callers
  that previously enabled `test-utils` should add
  `lambert_izzo_test_support` as a dev-dependency and import from
  `lambert_izzo_test_support::kepler::propagate` directly.
- **`pub mod test_utils`** and `src/test_helpers.rs` removed alongside
  the feature; the redundant external `tests/test_utils_smoke.rs` is
  gone (its coverage is fully subsumed by `src/lib.rs`'s inline tests).

### Added

- **`lambert_izzo_test_support` crate** (workspace-internal) consolidates
  dev fixtures previously duplicated across examples, benches, and
  integration tests:
  - `bodies::{MU_EARTH, MU_SUN, AU}` — standard astrodynamics constants
  - `rand_unit_vec(rng)` — rejection-sampling unit-vector helper
  - `kepler::propagate(r, v, dt, mu)` — universal-variable propagator
- **Doc-rubric coverage** on `solve_with_diagnostics`, `lambert_both_ways`,
  and the WASM wrappers (`solve_lambert`, `solve_lambert_request`).
  Each now has explicit `# Invariants` and `# Validity` sections that
  defer to `lambert`'s rubric.

### Changed (internal)

- **`vec3.rs` trimmed** from 9 functions to 4 (`dot`, `cross`, `norm`,
  `scale`). `add`, `sub`, `normalize`, `try_normalize`, `norm_squared`
  are inlined at their few remaining call sites in `geometry.rs` and
  the velocity reconstructor in `lib.rs`. The lib's inline tests now
  use `glam::DVec3` for symmetry with the examples and benches.
- **Hot-path: hoisted `y` out of the Householder / Halley iteration.**
  `tof::x_to_tof_with_y` and `tof::tof_derivatives_with_y` accept
  precomputed `y = sqrt(1 − λ²(1 − x²))`; the iteration loops compute
  it once per step and thread it into both calls, saving one
  `sqrt + arithmetic` per iteration on the dominant hot path.
- **Redundant `(1 − x²)` recomputation eliminated** in
  `tof::compute_psi` — now takes the precomputed value from its sole
  caller, `x_to_tof_lancaster`.
- **Examples, benches, and integration tests** dropped their
  hand-rolled vec helpers, local `MU_EARTH` constants (9 sites), and
  local `rand_unit_vec` impls (5 sites). All pull from
  `lambert_izzo_test_support` + `glam::DVec3`.
- **`glam`** promoted from `lambert_izzo`'s inline dev-dep to
  `[workspace.dependencies]`.

Behavior is preserved exactly — `examples/stress.rs` shows identical
solver-iteration histograms and conservation residuals at f64
round-off (~1e-12 single-rev, ~1e-14 multi-rev).

## [0.5.0] — 2026-04-26

### Changed (breaking)

- **Drop unit suffixes from public API.** Parameter and field names are
  now plain (`r1`, `r2`, `tof`, `mu`, `v1`, `v2`, `norm`) — the unit
  convention (km, s, km/s, km³/s²) lives in the docs, matching the
  prevailing Rust idiom. The crate is now dimensionally homogeneous in
  any consistent unit system; the SI choice is just a documentation
  default.
- **`LambertError::DegeneratePositionVector`** swaps its `which: u8`
  field for `position: Position`, where `Position` is a typed
  `R1 | R2` enum. Lets callers pattern-match on the variant instead of
  comparing magic numbers.
- **`NonFiniteParameter`** variants drop the `Km` suffix (`R1KmX` →
  `R1X`, etc.). `as_str()` returns `"r1.x"` etc.
- **`MIN_POSITION_NORM_KM`** renamed to `MIN_POSITION_NORM`.
- **`LambertError::NonPositiveTimeOfFlight::tof_s`** field renamed to
  `tof`; **`LambertError::NonPositiveMu::mu_km3_s2`** to `mu`;
  **`LambertError::DegeneratePositionVector::norm_km`** to `norm`.

### WASM adapter (`lambert_izzo_wasm` v0.4.0)

- Mirror types and request/response fields renamed to match the core
  crate. JS/TS callers now see `r1`, `r2`, `tof`, `mu`, `v1`, `v2`
  (camelCased per `serde(rename_all = "camelCase")`) instead of `r1Km`,
  `tofS`, etc.
- Adds a `PositionOutput` mirror enum for the core's new `Position`.

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
- Householder iteration with derivative-matched starters (Izzo Eq. 30, 31) and Halley `T_min` search for multi-rev infeasibility.
- Round-trip Kepler validation via universal-variable propagator;
  `examples/stress.rs` reproduces Izzo §5 statistical sweeps.
- Strict lint baseline: `clippy::pedantic` + bans on `unwrap`/`expect`/
  `panic`/`unreachable` in lib code.

[Unreleased]: https://github.com/sakobu/izzos_lambert/compare/v0.5.0...HEAD
[0.5.0]: https://github.com/sakobu/izzos_lambert/compare/v0.4.0...v0.5.0
[0.4.0]: https://github.com/sakobu/izzos_lambert/compare/v0.3.0...v0.4.0
[0.3.0]: https://github.com/sakobu/izzos_lambert/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/sakobu/izzos_lambert/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/sakobu/izzos_lambert/releases/tag/v0.1.0
