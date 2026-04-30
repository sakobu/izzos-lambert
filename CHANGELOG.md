# Changelog

All notable changes to `lambert_izzo` are documented here. The format
follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and
the project adheres to [Semantic Versioning](https://semver.org/) once
it reaches `1.0`.

## [Unreleased]

### Optional-feature examples

- **New `examples/batch.rs`** (gated on `rayon`) — drives `lambert_par`
  over 10 000 randomized Earth-scale inputs and reports wall-clock
  throughput plus the mean Householder iteration count across successful
  solves. Run with `cargo run --release --example batch --features rayon`.
- **New `examples/serde.rs`** (gated on `serde`) — round-trips a
  `LambertSolutions` and a `LambertError` through `serde_json`,
  asserting `PartialEq` equality on both ends. Run with
  `cargo run --release --example serde --features serde`.

### Type-enforced revolution cap (breaking)

The `MAX_MULTI_REV_PAIRS = 32` cap on multi-rev pairs is now type-level
rather than enforced by silent runtime clamping.

- **New `BoundedRevs` type** — newtype around `NonZeroU32` constrained to
  `1..=BoundedRevs::MAX`. Constructed via `BoundedRevs::try_new(u32) ->
  Result<Self, RevsOutOfRange>`. The MAX constant is statically asserted
  to match `MAX_MULTI_REV_PAIRS`.
- **`RevolutionBudget::UpTo` now wraps `BoundedRevs`** (was
  `NonZeroU32`).
- **`RevolutionBudget::up_to`** changed signature: takes
  `BoundedRevs` (total), no longer takes `u32` and silently collapses
  `0` to `SingleOnly`. Passing `0` was always a footgun.
- **New `RevolutionBudget::try_up_to(u32) -> Result<Self,
  RevsOutOfRange>`** — ergonomic fallible constructor for the common
  case (literal or external `u32`). Use `RevolutionBudget::SingleOnly`
  directly to skip multi-rev.
- **`RevolutionBudget::max()` returns `Option<BoundedRevs>`** (was
  `u32`). `None` for `SingleOnly`, `Some(b)` for `UpTo(b)`. Removes the
  ambiguous `0` return and matches the type-honest direction of
  `BoundedRevs`.
- **New `RevolutionBudget::iter_revs() -> impl Iterator<Item = BoundedRevs>`** —
  yields validated `BoundedRevs` values `1..=b` for `UpTo(b)` and is
  empty for `SingleOnly`. Canonical way to drive a multi-rev loop
  without re-materializing the `0` upper-bound sentinel; emitting
  `BoundedRevs` rather than raw `u32` keeps the `1..=BoundedRevs::MAX`
  invariant inside the type system all the way to the
  `MultiRevPair::n_revs` field on the returned solution.
- **`MultiRevPair::n_revs` and `MultiRevPairDiagnostics::n_revs` are now
  `BoundedRevs`** (was `u32`); callers needing the raw count use
  `.n_revs.get()`. Closes the last `u32`-shaped invariant on the
  multi-rev public API. JS callers are unaffected — the WASM adapter
  still surfaces `nRevs: number`.
- **`BoundedRevs` gains `Display`, `PartialOrd`/`Ord`, and `Hash`
  derives** — pure additions. Lets `println!("M={}", pair.n_revs)`,
  ordering checks (`pair.n_revs > prev`), and use as a hash-map key
  drop in without `.get()` round-trips.
- **New `RevsOutOfRange` error type** — standalone, `thiserror::Error`,
  serde-compatible. Distinct from `LambertError` because it represents
  construction-time validation, not solver-runtime failure.
- **Kernel simplified** — the `m_max.min(MAX_MULTI_REV_PAIRS)` clamp in
  `find_xy` is gone; the type now carries the invariant.

### Multi-rev silent-skip diagnostic

- **New `LambertSolutions::max_feasible_revs()` method** — returns
  `Option<BoundedRevs>` carrying the highest revolution count `M` for
  which a feasible `(long_period, short_period)` pair was found at the
  requested TOF. `None` when `RevolutionBudget::SingleOnly` was used or
  no multi-rev branch was feasible. Lets callers programmatically detect
  the silent-skip behavior at the `T_min(M) > tof` boundary that until
  now was only documented in `lambert`'s validity rubric. Pairs
  symmetrically with `RevolutionBudget::max()` — both expose "highest
  relevant `M`" through the same type. JS callers derive the same value
  from `response.multi.at(-1)?.nRevs ?? null`.

### WASM adapter (breaking)

- **`LambertRequest.maxRevs`** changes type from `number` to
  `number | null`. Pass `null` to skip multi-rev (formerly `0`); pass
  `1..=32` to search multi-rev branches. Out-of-range values reject
  with a new `LambertErrorOutput::RevsOutOfRange { requested, max }`
  variant.
- **New `From<RevsOutOfRange> for LambertErrorOutput`** impl wires the
  core's validation error into the JS-facing tagged union.

## [1.0.0] — 2026-04-30

The first stable release. Workspace shares `version = "1.0.0"` across all
three crates; subsequent releases bump them together.

### Public API consolidated to two entry points (breaking)

The core crate's six entry points collapse to two:

- **`lambert(&LambertInput) -> Result<LambertSolutions, LambertError>`** —
  canonical single solve. Takes a `&LambertInput` instead of six
  positional arguments.
- **`lambert_par(&[LambertInput])`** (`rayon` feature) — parallel batch
  solve, returns an `IndexedParallelIterator` over per-input results.

Removed:

- `solve_with_diagnostics(...)` — diagnostics now live inside every
  `LambertSolutions` (see below), so the dichotomy is gone.
- `lambert_both_ways(...)` — saved no work (geometry was recomputed
  internally for short vs. long); call `lambert` twice with the two
  `TransferWay` values.
- `lambert_iter(&[LambertInput])` — was a one-line `.iter().map(...)`.
- `LambertInput::solve()` — the free function `lambert(&input)` is the
  single canonical entry point.
- `BothWaysSolutions` (the return type of `lambert_both_ways`).

Renamed:

- `lambert_par_iter` → `lambert_par`.

### Diagnostics live inside `LambertSolutions` (breaking)

Every `lambert(...)` call returns `LambertSolutions { single, multi,
diagnostics }`. The new `diagnostics` field carries the per-branch
[`SolverDiagnostics`] (Householder iteration count) for the `single`
branch and every multi-rev pair.

`SolverDiagnostics::lancaster_blanchard_x` was removed — it exposed a
kernel-internal variable that carried no consumer signal (a converged
`Ok(...)` already implies tolerance was met, and `LambertError::NoConvergence`
already carries the failure information).

### Kernel internals hidden (breaking)

- `pub mod constants` is now private. The Householder / Halley
  tolerances and iteration caps are not part of the public API contract.
- `MAX_MULTI_REV_PAIRS` (the only externally relevant constant — it
  bounds the public return type) stays at the crate root.

### `MultiRevSet` / `MultiRevDiagnostics` newtypes (breaking)

`LambertSolutions::multi` and `LambertDiagnostics::multi` are now newtypes
wrapping the underlying `ArrayVec`. Consumers see only
`Deref<Target = [MultiRevPair]>` (or `[MultiRevPairDiagnostics]`) and
`IntoIterator`. The `arrayvec` crate no longer leaks into the public API.

JSON serialization (under the `serde` feature) uses `#[serde(transparent)]`
so the on-the-wire shape stays a flat array.

### Tests reorganized

`#[cfg(test)] mod tests` moved out of `lib.rs` into per-scenario
submodules under `src/tests/`: `single_rev`, `multi_rev`, `errors`,
`regimes`, `interop`, `kepler_roundtrip`. Test bodies are unchanged
modulo the API consolidation above.

### WASM adapter (`lambert_izzo_wasm` v1.0.0, breaking)

- `LambertSolutionOutput.diagnostics` field removed.
- `LambertResponse.diagnostics: LambertDiagnosticsOutput` added (mirrors
  the core's top-level `diagnostics`).
- `SolverDiagnosticsOutput.x` field removed.
- New `solveLambertBatch(requests)` that returns `Vec<BatchResult>` —
  per-input tagged `{ kind: "ok", response } | { kind: "err", error }`.
- New `MultiRevPairDiagnosticsOutput` and `LambertDiagnosticsOutput`
  mirror types.
- `LambertErrorOutput::Unknown { message }` retained as the forward-compat
  fallback (the core's `LambertError` stays `#[non_exhaustive]`).
- Real `README.md` and a single-page browser demo at `examples/web/`.

### Workspace versioning

- `[workspace.package]` block declares `version = "1.0.0"`; all three
  member crates use `version.workspace = true`.
- `lambert_izzo_test_support` continues as `publish = false`
  (workspace-internal).

## [0.5.0] — 2026-04-26 (pre-1.0)

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
  and the WASM wrappers. Each now has explicit `# Invariants` and
  `# Validity` sections that defer to `lambert`'s rubric.

### Changed (internal)

- **`vec3.rs` trimmed** from 9 functions to 4 (`dot`, `cross`, `norm`,
  `scale`). `add`, `sub`, `normalize`, `try_normalize`, `norm_squared`
  are inlined at their few remaining call sites.
- **Hot-path: hoisted `y` out of the Householder / Halley iteration.**
- **Redundant `(1 − x²)` recomputation eliminated** in `tof::compute_psi`.
- **Examples, benches, and integration tests** dropped their hand-rolled
  vec helpers, local `MU_EARTH` constants (9 sites), and local
  `rand_unit_vec` impls (5 sites). All pull from `lambert_izzo_test_support`
  + `glam::DVec3`.
- **`glam`** promoted from `lambert_izzo`'s inline dev-dep to
  `[workspace.dependencies]`.

### Public-API renames (breaking, pre-1.0)

- **Drop unit suffixes from public API.** Parameter and field names are
  now plain (`r1`, `r2`, `tof`, `mu`, `v1`, `v2`, `norm`).
- **`LambertError::DegeneratePositionVector`** swaps its `which: u8`
  field for `position: Position`, where `Position` is a typed
  `R1 | R2` enum.
- **`NonFiniteParameter`** variants drop the `Km` suffix (`R1KmX` → `R1X`,
  etc.).
- **`MIN_POSITION_NORM_KM`** renamed to `MIN_POSITION_NORM`.
- **`LambertError::NonPositiveTimeOfFlight::tof_s`** field renamed to
  `tof`; **`NonPositiveMu::mu_km3_s2`** to `mu`;
  **`DegeneratePositionVector::norm_km`** to `norm`.

### WASM adapter (`lambert_izzo_wasm` v0.4.0)

- Mirror types and request/response fields renamed to match the core
  crate. JS/TS callers now see `r1`, `r2`, `tof`, `mu`, `v1`, `v2`
  (camelCased per `serde(rename_all = "camelCase")`) instead of `r1Km`,
  `tofS`, etc.
- Adds a `PositionOutput` mirror enum for the core's new `Position`.

## [0.4.0] — 2026-04-26 (pre-1.0)

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

## [0.3.0] — 2026-04-26 (pre-1.0)

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

## [0.2.0] — 2026-04-26 (pre-1.0)

- Workspace restructure: `crates/lambert_izzo` (core) and
  `crates/lambert_izzo_wasm` (`wasm-bindgen` adapter).
- Documentation polish on `pub fn lambert`; tightened `Cargo.toml`
  metadata for `crates.io` publication.

## [0.1.0] — 2026-04-26 (pre-1.0)

- Initial implementation of Izzo's revisited Lambert solver
  (single + multi-revolution, short + long way).
- Three-regime time-of-flight dispatch (Battin / Lancaster–Blanchard /
  Lagrange) with documented thresholds in `constants.rs`.
- Householder iteration with derivative-matched starters (Izzo Eq. 30, 31) and Halley `T_min` search for multi-rev infeasibility.
- Round-trip Kepler validation via universal-variable propagator;
  `examples/stress.rs` reproduces Izzo §5 statistical sweeps.
- Strict lint baseline: `clippy::pedantic` + bans on `unwrap`/`expect`/
  `panic`/`unreachable` in lib code.

[Unreleased]: https://github.com/sakobu/izzos_lambert/compare/v1.0.0...HEAD
[1.0.0]: https://github.com/sakobu/izzos_lambert/compare/v0.5.0...v1.0.0
[0.5.0]: https://github.com/sakobu/izzos_lambert/compare/v0.4.0...v0.5.0
[0.4.0]: https://github.com/sakobu/izzos_lambert/compare/v0.3.0...v0.4.0
[0.3.0]: https://github.com/sakobu/izzos_lambert/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/sakobu/izzos_lambert/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/sakobu/izzos_lambert/releases/tag/v0.1.0
