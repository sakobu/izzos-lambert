# lambert_izzo

A Rust port of Dario Izzo's revisited Lambert solver from the 2014 paper
_"Revisiting Lambert's Problem"_
([arXiv:1403.2705](https://arxiv.org/abs/1403.2705) / Celestial Mechanics &
Dynamical Astronomy). A local copy of the paper lives at
[`docs/izzo.pdf`](docs/izzo.pdf); for a friendly intro to the problem,
read [`docs/concepts.md`](docs/concepts.md) first.

> Pre-1.0 — the API and feature surface are stable enough for production
> use but may shift before `1.0`. Breaking changes are tracked in
> [`CHANGELOG.md`](CHANGELOG.md).

## What this solves (and what it doesn't)

| Capability                                                             | Supported |
| ---------------------------------------------------------------------- | --------- |
| Single-revolution two-body transfers (elliptic, parabolic, hyperbolic) | ✅        |
| Multi-revolution transfers, both long- and short-period branches       | ✅        |
| Short-way and long-way arc selection                                   | ✅        |
| Frame-invariant inputs/outputs (any inertial frame the caller chooses) | ✅        |
| `no_std` and WASM (`wasm32-unknown-unknown`) builds                    | ✅        |
| Optional Rayon-backed parallel batch evaluation                        | ✅        |
| Patched-conic / sphere-of-influence transitions                        | ❌ caller |
| Lunar swing-by, n-body perturbations, J2/J3 effects                    | ❌ caller |
| Low-thrust / continuous-thrust transfers                               | ❌        |
| Outer-loop optimization (porkchop-min, primer-vector, …)               | ❌ caller |

"❌ caller" means this is the right Lambert solver to call from inside
those higher-level routines — but the routine itself isn't part of this
crate.

Supports:

- Single-revolution transfers
- Multi-revolution transfers (long-period and short-period branches)
- Short-way and long-way transfers (`TransferWay::Short` /
  `TransferWay::Long`, or `lambert_both_ways(...)` for one call returning
  both); prograde vs retrograde is the caller's choice via the `(r1, r2)`
  ordering, since `r1 × r2` defines the resulting orbit's
  angular-momentum direction
- Hyperbolic transfers on the single-rev branch
- `no_std`-friendly — pulls only `arrayvec`, `num-traits` (with `libm`),
  and `thiserror` (`std`-feature off) at runtime
- WASM-compatible math kernel
  (`cargo build --target wasm32-unknown-unknown --no-default-features --lib`)
- Zero hard math-library dependency — public surface is `[f64; 3]`

## Features

| Feature | Default | Effect                                                                                                              |
| ------- | ------- | ------------------------------------------------------------------------------------------------------------------- |
| `serde` | off     | Adds `Serialize`/`Deserialize` derives on every public type, including `LambertError`.                              |
| `rayon` | off     | Enables `lambert_par_iter` for parallel batch evaluation. Pulls in `std` transitively — incompatible with `no_std`. |

> Need a Kepler propagator to round-trip Lambert solutions in your own
> tests? It used to live behind a `test-utils` feature here; it now sits
> in the workspace-internal `lambert_izzo_test_support` crate
> (`publish = false`). External consumers should vendor or reimplement
> the propagator — it's ~80 lines of universal-variable / Stumpff math.

MSRV: **Rust 1.85** (the first release with edition 2024 stable).

## Units

The crate is unit-agnostic at the type level — every quantity is plain
`f64` or `[f64; 3]`. The convention used in the docs and examples is SI
for astrodynamics work:

| Quantity                | Unit   |
| ----------------------- | ------ |
| Position                | km     |
| Velocity                | km/s   |
| Time of flight          | s      |
| Gravitational parameter | km³/s² |

Any consistent unit system works — pass `r` in meters and `mu` in m³/s²
and you get velocities in m/s. The math is dimensionally homogeneous;
unit safety is the caller's responsibility, helped by your own variable
names (`r1_eci_km`, `mu_earth_km3_s2`, etc.) at the call site.

The algorithm is also frame-invariant under any inertial frame — pass
`r1`, `r2` in the same inertial frame (ECI, HCRS, MCI, …) and the
returned velocities are in that same frame. The function signature is
frame-agnostic; the calling code's variable names carry the frame info.

## Usage

```rust
use lambert_izzo::{lambert, RevolutionBudget, TransferWay};

// LEO → MEO Hohmann transfer.
let mu = 398_600.441_8;
let r1 = [7000.0, 0.0, 0.0];
let r2 = [-12_000.0, 1.0, 0.0]; // 1 km off-axis avoids colinearity
let a = (7000.0 + 12_000.0) / 2.0;
let tof = std::f64::consts::PI * (a.powi(3) / mu).sqrt();

let solutions = lambert(
    r1, r2, tof, mu,
    TransferWay::Short, RevolutionBudget::SingleOnly,
).unwrap();
let v1 = solutions.single.v1;
let v2 = solutions.single.v2;
```

The signature:

```rust
pub fn lambert(
    r1: [f64; 3],                  // initial position, any inertial frame
    r2: [f64; 3],                  // final position, same frame
    tof: f64,                      // time of flight, > 0
    mu: f64,                       // gravitational parameter, > 0
    way: TransferWay,              // Short or Long way around the transfer plane
    revolutions: RevolutionBudget, // SingleOnly or UpTo(NonZero<u32>)
) -> Result<LambertSolutions, LambertError>
```

The returned `LambertSolutions` is a typed shape — single-revolution
always present, multi-revolution branches in `(long_period, short_period)`
pairs:

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
    pub v1: [f64; 3],
    pub v2: [f64; 3],
}
```

For the iteration count and Lancaster–Blanchard `x` (useful for
diagnosing multi-rev branches), use `solve_with_diagnostics`:

```rust
use lambert_izzo::solve_with_diagnostics;
let (sols, diag) = solve_with_diagnostics(
    r1, r2, tof, mu,
    TransferWay::Short, RevolutionBudget::up_to(3),
)?;
println!("converged in {} iters", diag.single.iters);
```

For the porkchop-plot pattern (you want both ways), use
`lambert_both_ways`:

```rust
use lambert_izzo::lambert_both_ways;
let both = lambert_both_ways(r1, r2, tof, mu, RevolutionBudget::up_to(3))?;
let short_v1 = both.short.single.v1;
let long_v1 = both.long.single.v1;
```

`LambertError` is a `thiserror` enum with structured fields:

```rust
match lambert(r1, r2, tof, mu, TransferWay::Short, RevolutionBudget::SingleOnly) {
    Ok(sols) => /* … */,
    Err(LambertError::NonFiniteInput { parameter, value })             => /* … */,
    Err(LambertError::NonPositiveTimeOfFlight { tof })                  => /* … */,
    Err(LambertError::NonPositiveMu { mu })                             => /* … */,
    Err(LambertError::DegeneratePositionVector { position, norm })      => /* … */,
    Err(LambertError::CollinearGeometry { sin_angle })                  => /* … */,
    Err(LambertError::NoConvergence { iterations, last_step, n_revs })  => /* … */,
    Err(_) => /* … */,
}
```

### Math-library interop

The crate has no hard math-library dependency. Both
`nalgebra::Vector3<f64>` and `glam::DVec3` already convert to/from
`[f64; 3]` natively, so callers using either library can pass and
receive vectors without any feature flag:

```rust
// nalgebra:
let r1: [f64; 3] = nalgebra::Vector3::new(7000.0, 0.0, 0.0).into();
let v1: nalgebra::Vector3<f64> = solutions.single.v1.into();

// glam:
let r2 = glam::DVec3::new(0.0, 7000.0, 0.0).to_array();
let v2 = glam::DVec3::from_array(solutions.single.v2);
```

## Validation

The `stress` example runs 100,000 random Earth-orbit geometries each for
single-rev and multi-rev (up to `M = 5`), checking vis-viva energy and
angular-momentum conservation between the returned `(v1, v2)` pair.
Random ranges:

- Single-rev: `r ∈ [3500, 28_000]` km, `tof ∈ [100, 50_000]` s
- Multi-rev: `r ∈ [5600, 10_500]` km, `tof ∈ [10_000, 250_000]` s

| Mode       | Convergence | Avg iters | Paper avg | Max iters | Max ΔE/E | Max ΔL/L |
| ---------- | ----------- | --------- | --------- | --------- | -------- | -------- |
| Single-rev | 100%        | 2.083     | 2.1       | 3         | 1.44e-12 | 4.14e-12 |
| Multi-rev  | 100%        | 2.992     | 3.3       | 7         | 3.02e-14 | 5.63e-14 |

Iteration counts match the paper's reported figures; conservation
residuals sit at f64 round-off.

## Performance

Criterion benchmarks under `crates/lambert_izzo/benches/`. Numbers below
are from an Apple Silicon laptop (release profile, single thread except
for the parallel batch row).

| Workload                                      | Throughput     | Per call |
| --------------------------------------------- | -------------- | -------- |
| Single-rev (random Earth orbits)              | ~2.7 M calls/s | ~360 ns  |
| Multi-rev `M=3` (Earth orbits)                | ~735 K calls/s | ~1.4 µs  |
| Sequential batch via `lambert_iter`           | ~1.1 M calls/s | ~900 ns  |
| Parallel batch via `lambert_par_iter` (rayon) | ~8.6 M calls/s | ~116 ns  |

The parallel batch shows ~7.7× speedup over sequential on this machine.
Reproduce with:

```
cargo bench --bench single_rev
cargo bench --bench multi_rev
cargo bench --bench batch --features rayon
```

## Batch / streaming API

For porkchop plots, multi-shooter loops, or any workload with thousands
of Lambert calls, `lambert_iter` and (under the `rayon` feature)
`lambert_par_iter` map a slice of `LambertInput`:

```rust
use lambert_izzo::{LambertInput, RevolutionBudget, TransferWay, lambert_iter};

let inputs: Vec<LambertInput> = (0..10_000)
    .map(|_| LambertInput { /* … */ # r1: [7000.0, 0.0, 0.0],
        # r2: [0.0, 7000.0, 0.0], tof: 1500.0, mu: 398_600.4418,
        # way: TransferWay::Short,
        # revolutions: RevolutionBudget::SingleOnly,
    })
    .collect();

let total_dv: f64 = lambert_iter(&inputs)
    .filter_map(Result::ok)
    .map(|sols| sols.single.v1[0].abs())
    .sum();
```

With `--features rayon`, swap `lambert_iter` for `lambert_par_iter` and
the work fans out across the thread pool.

## Building

```
cargo build --release
cargo test --release
cargo run --release --example demo
cargo run --release --example stress
cargo run --release --example errors
```

The `errors` example walks every `LambertError` variant — useful as a
template for caller-side error handling.

Toolchain pinned via `rust-toolchain.toml` (1.88.0) for development;
MSRV declared in `Cargo.toml` is 1.85. Edition 2024. Runtime
dependencies are [`thiserror`](https://docs.rs/thiserror) (no_std mode)
for the error type, [`arrayvec`](https://docs.rs/arrayvec) (no_std) for
the bounded multi-rev return, and
[`num-traits`](https://docs.rs/num-traits) (with `libm`) for `no_std`
math.

## Implementation notes

- Modular layout under `src/`:
  - `constants.rs` — every named tolerance with rationale.
  - `error.rs` — structured `LambertError` enum.
  - `vec3.rs` — minimal internal `[f64; 3]` helpers: `dot`, `cross`,
    `norm`, `scale`. Trivial component-wise operations are inlined at
    their call sites.
  - `geometry.rs` — chord, semi-perimeter, λ, transfer-plane basis;
    constructed once per call.
  - `tof.rs` — three-regime TOF dispatch + analytic derivatives (Eq. 22).
    The `_with_y` variants accept a precomputed `y = √(1 − λ²(1 − x²))`
    so the Householder loop computes it once per step instead of twice.
  - `root_finding.rs` — Householder (Eq. 30/31 starters) + Halley
    `T_min` search.
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
