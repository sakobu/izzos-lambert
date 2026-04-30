# Changelog

All notable changes to `lambert_izzo` are documented here. The format
follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and
the project follows [Semantic Versioning](https://semver.org/).

## [Unreleased]

## [1.0.0] — 2026-04-30

Initial public release. `lambert_izzo` is a `no_std`-friendly Rust
implementation of Izzo's revisited Lambert solver
([arXiv:1403.2705](https://arxiv.org/abs/1403.2705)).

### Solver

- Single-revolution and multi-revolution (short- and long-way) Lambert
  solves via `lambert(&LambertInput) -> Result<LambertSolutions, LambertError>`.
- Three-regime time-of-flight dispatch — Battin hypergeometric series
  near the parabolic point, Lancaster–Blanchard in the transitional
  zone, Lagrange elsewhere — with all thresholds documented in
  `constants.rs`.
- Cubic Householder iteration with derivative-matched starters
  (Izzo Eq. 30, 31) and Halley `T_min` search for multi-revolution
  feasibility detection.
- Round-trip Kepler validation via a universal-variable Stumpff
  propagator; `examples/stress.rs` runs 200,000 randomized trials
  with vis-viva and angular-momentum conservation checks.

### Public API

- `[f64; 3]` for all positions and velocities — no hard math-library
  dependency. `nalgebra::Vector3<f64>` and `glam::DVec3` convert
  natively.
- Unit-agnostic: pass any consistent unit system. Documented SI
  convention is km / km/s / s / km³/s².
- Frame-invariant: the solver returns velocities in whatever inertial
  frame the inputs are in (ECI, HCRS, MCI, ...). Frame info lives at
  the call site, not in field names.
- `LambertError` is a flat `thiserror` enum with structured fields
  (`Position`, `NonFiniteParameter`, `RevsOutOfRange`); zero
  `unwrap`/`expect`/`panic!` in lib code, mechanically enforced by
  clippy denials.
- Multi-rev capacity is bounded at compile time via `BoundedRevs`
  (newtype around `NonZeroU32`, `1..=32`); the cap is statically
  asserted to match `MAX_MULTI_REV_PAIRS`. No heap allocation on the
  solver hot path.
- Diagnostics (Householder iteration counts) ride along on every
  `LambertSolutions` — no separate "with-diagnostics" entry point.

### Cargo features (off by default)

- **`serde`** — derives `Serialize` / `Deserialize` on every public
  type, including `LambertError`. Stays `no_std`-compatible.
- **`rayon`** — enables `lambert_par(&[LambertInput])`, a parallel
  batch solver. Pulls `std` transitively; not `no_std`-compatible.

### Workspace layout

- **`lambert_izzo`** — the canonical Rust crate. `no_std`-friendly,
  `[f64; 3]`-only public surface.
- **`lambert_izzo_wasm`** — `wasm-bindgen` adapter exposing
  JavaScript / TypeScript request and response types via `tsify`.
  Errors map to a tag-discriminated union; a forward-compat
  `Unknown { message }` variant guards against future core variants.
- **`lambert_izzo_test_support`** — workspace-internal dev fixtures
  (`publish = false`): SI body constants, rejection-sampling unit
  vectors, deterministic batch generator, universal-variable Kepler
  propagator. Used as a `[dev-dependencies]` path crate by the core's
  examples, benches, and tests.

### Tooling

- MSRV: Rust 1.85 (first edition-2024 stable).
- Examples: `demo`, `stress`, `errors`, `batch` (`rayon`-gated), and
  `serde` (`serde`-gated).
- Benches: `single_rev`, `multi_rev`, `batch` via Criterion.
- Dual-licensed `MIT OR Apache-2.0`.
