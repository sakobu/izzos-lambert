# lambert_izzo

A Rust port of Dario Izzo's revisited Lambert solver from the 2014 paper
_"Revisiting Lambert's Problem"_
([arXiv:1403.2705](https://arxiv.org/abs/1403.2705) / Celestial Mechanics &
Dynamical Astronomy). A local copy of the paper lives at
[`docs/izzo.pdf`](docs/izzo.pdf).

Supports:

- Single-revolution transfers
- Multi-revolution transfers (long-period and short-period branches)
- Short-way and long-way transfers (`TransferWay::Short` / `TransferWay::Long`,
  or `lambert_both_ways(...)` for one call returning both); prograde vs retrograde
  is the caller's choice via the `(r1, r2)` ordering, since
  `r1 × r2` defines the resulting orbit's angular-momentum direction
- Hyperbolic transfers on the single-rev branch
- `no_std`-friendly — pulls only `arrayvec`, `num-traits` (with `libm`), and
  `thiserror` (`std`-feature off) at runtime
- WASM-compatible math kernel (`cargo build --target wasm32-unknown-unknown --no-default-features --lib`)
- Zero hard math-library dependency — public surface is `[f64; 3]`

## Features

| Feature      | Default | Effect                                                                                                |
| ------------ | ------- | ----------------------------------------------------------------------------------------------------- |
| `serde`      | off     | Adds `Serialize`/`Deserialize` derives on every public type, including `LambertError`.               |
| `test-utils` | off     | Promotes the universal-variable Kepler propagator to `lambert_izzo::test_utils::kepler_propagate` so downstream integration tests can round-trip-validate Lambert solutions without re-implementing it. |

MSRV: **Rust 1.85** (the first release with edition 2024 stable).

## Units

Public API uses SI conventions:

| Quantity                | Suffix    | Unit   |
| ----------------------- | --------- | ------ |
| Position                | `_km`     | km     |
| Velocity                | `_km_s`   | km/s   |
| Time                    | `_s`      | s      |
| Gravitational parameter | `_km3_s2` | km³/s² |

The algorithm is mathematically frame-invariant under any inertial frame —
pass `r1_km`, `r2_km` in the same inertial frame (ECI, HCRS, MCI, …) and
the returned velocities are in that same frame. The function signature
itself is frame-agnostic; the calling code's variable names carry the
frame info.

## Usage

```rust
use lambert_izzo::{lambert, RevolutionBudget, TransferWay};

// LEO → MEO Hohmann transfer.
let mu_km3_s2 = 398_600.441_8;
let r1_km = [7000.0, 0.0, 0.0];
let r2_km = [-12_000.0, 1.0, 0.0]; // 1 km off-axis avoids colinearity
let a_km = (7000.0 + 12_000.0) / 2.0;
let tof_s = std::f64::consts::PI * (a_km.powi(3) / mu_km3_s2).sqrt();

let solutions = lambert(
    r1_km, r2_km, tof_s, mu_km3_s2,
    TransferWay::Short, RevolutionBudget::SingleOnly,
).unwrap();
let v1_km_s = solutions.single.v1_km_s;
let v2_km_s = solutions.single.v2_km_s;
```

The signature:

```rust
pub fn lambert(
    r1_km: [f64; 3],        // initial position (km), any inertial frame
    r2_km: [f64; 3],        // final position (km), same frame
    tof_s: f64,             // time of flight (s), > 0
    mu_km3_s2: f64,         // gravitational parameter (km³/s²), > 0
    way: TransferWay,       // Short or Long way around the transfer plane
    revolutions: RevolutionBudget, // SingleOnly or UpTo(NonZero<u32>)
) -> Result<LambertSolutions, LambertError>
```

The returned `LambertSolutions` is a typed shape — single-revolution always
present, multi-revolution branches in `(long_period, short_period)` pairs:

```rust
pub struct LambertSolutions {
    pub single: LambertSolution,
    pub multi: ArrayVec<MultiRevPair, MAX_MULTI_REV_PAIRS>, // up to 32, stack-allocated
}

pub struct MultiRevPair {
    pub n_revs: u32,
    pub long_period: LambertSolution,
    pub short_period: LambertSolution,
}

pub struct LambertSolution {
    pub v1_km_s: [f64; 3],
    pub v2_km_s: [f64; 3],
}
```

For the iteration count and Lancaster–Blanchard `x` (useful for diagnosing
multi-rev branches), use `solve_with_diagnostics`:

```rust
use lambert_izzo::solve_with_diagnostics;
let (sols, diag) = solve_with_diagnostics(
    r1_km, r2_km, tof_s, mu_km3_s2,
    TransferWay::Short, RevolutionBudget::up_to(3),
)?;
println!("converged in {} iters", diag.single.iters);
```

For the porkchop-plot pattern (you want both ways), use `lambert_both_ways`:

```rust
use lambert_izzo::lambert_both_ways;
let both = lambert_both_ways(r1_km, r2_km, tof_s, mu_km3_s2, RevolutionBudget::up_to(3))?;
let short_v1 = both.short.single.v1_km_s;
let long_v1 = both.long.single.v1_km_s;
```

`LambertError` is a `thiserror` enum with structured fields:

```rust
match lambert(r1_km, r2_km, tof_s, mu_km3_s2, TransferWay::Short, RevolutionBudget::SingleOnly) {
    Ok(sols) => /* … */,
    Err(LambertError::NonFiniteInput { parameter, value }) => /* … */,
    Err(LambertError::NonPositiveTimeOfFlight { tof_s })       => /* … */,
    Err(LambertError::NonPositiveMu { mu_km3_s2 })             => /* … */,
    Err(LambertError::DegeneratePositionVector { which, norm_km }) => /* … */,
    Err(LambertError::CollinearGeometry { sin_angle })         => /* … */,
    Err(LambertError::NoConvergence { iterations, last_step, n_revs }) => /* … */,
    Err(_) => /* … */,
}
```

### Math-library interop

The crate has no hard math-library dependency. Both `nalgebra::Vector3<f64>`
and `glam::DVec3` already convert to/from `[f64; 3]` natively, so callers
using either library can pass and receive vectors without any feature flag:

```rust
// nalgebra:
let r1_km: [f64; 3] = nalgebra::Vector3::new(7000.0, 0.0, 0.0).into();
let v1_na: nalgebra::Vector3<f64> = solutions.single.v1_km_s.into();

// glam:
let r2_km = glam::DVec3::new(0.0, 7000.0, 0.0).to_array();
let v2_glam = glam::DVec3::from_array(solutions.single.v2_km_s);
```

## Validation

The `stress` example runs 100,000 random Earth-orbit geometries each for
single-rev and multi-rev (up to `M = 5`), checking vis-viva energy and
angular-momentum conservation between the returned `(v1_km_s, v2_km_s)`
pair. Random ranges:

- Single-rev: `r ∈ [3500, 28_000]` km, `tof ∈ [100, 50_000]` s
- Multi-rev: `r ∈ [5600, 10_500]` km, `tof ∈ [10_000, 250_000]` s

| Mode       | Convergence | Avg iters | Paper avg | Max iters | Max ΔE/E | Max ΔL/L |
| ---------- | ----------- | --------- | --------- | --------- | -------- | -------- |
| Single-rev | 100%        | 2.084     | 2.1       | 3         | 8.41e-12 | 1.86e-12 |
| Multi-rev  | 100%        | 2.992     | 3.3       | 7         | 3.00e-14 | 1.37e-13 |

Iteration counts match the paper's reported figures.

## Building

```
cargo build --release
cargo test --release
cargo run --release --example demo
cargo run --release --example stress
```

Toolchain pinned via `rust-toolchain.toml` (1.88.0) for development; MSRV
declared in `Cargo.toml` is 1.85. Edition 2024. Runtime dependencies are
[`thiserror`](https://docs.rs/thiserror) (no_std mode) for the error type,
[`arrayvec`](https://docs.rs/arrayvec) (no_std) for the bounded multi-rev
return, and [`num-traits`](https://docs.rs/num-traits) (with `libm`) for
`no_std` math.

## Implementation notes

- Modular layout under `src/`:
  - `constants.rs` — every named tolerance with rationale.
  - `error.rs` — structured `LambertError` enum.
  - `vec3.rs` — internal `[f64; 3]` helpers (dot, cross, norm, etc.).
  - `geometry.rs` — chord, semi-perimeter, λ, transfer-plane basis;
    constructed once per call.
  - `tof.rs` — three-regime TOF dispatch + analytic derivatives (Eq. 22).
  - `root_finding.rs` — Householder (Eq. 30/31 starters) + Halley `T_min`
    search.
  - `lib.rs` — public API + integration tests.
- TOF evaluation blends Battin's series (Eq. 20) for `|x − 1| ≤ 0.01`,
  Lancaster–Blanchard (Eq. 18) for `0.01 < |x − 1| ≤ 0.2`, and Lagrange
  (Eq. 9) elsewhere. The Battin path uses a direct series sum of
  `2F1(3, 1; 5/2; S1)`.
- Root finding uses Householder's 3rd-order method per the paper, with
  separate tolerances `1e-5` for `M = 0` and `1e-8` for `M > 0`.
- For multi-rev, `T_min` is found via Halley's method on `dT/dx = 0`
  before deciding which revolution counts admit solutions.
- Initial guesses follow Eq. 30 (single-rev) and Eq. 31 (multi-rev).
- Velocity reconstruction follows Algorithm 1.
- Multi-rev branches are clamped at `MAX_MULTI_REV_PAIRS = 32`. The
  bounded `ArrayVec` return means `lambert(...)` does no heap allocation
  on any code path (`MAX_MULTI_REV_PAIRS * sizeof(MultiRevPair)` lives
  on the stack).

## License

MIT OR Apache-2.0
