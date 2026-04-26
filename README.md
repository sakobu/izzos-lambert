# lambert_izzo

A Rust port of Dario Izzo's revisited Lambert solver from the 2014 paper
_"Revisiting Lambert's Problem"_
([arXiv:1403.2705](https://arxiv.org/abs/1403.2705) / Celestial Mechanics &
Dynamical Astronomy). A local copy of the paper lives at
[`docs/izzo.pdf`](docs/izzo.pdf).

Supports:

- Single-revolution transfers
- Multi-revolution transfers (long-period and short-period branches)
- Short-way and long-way transfers (`way: TransferWay`); prograde vs retrograde
  is the caller's choice via the `(r1, r2)` ordering, since
  `r1 × r2` defines the resulting orbit's angular-momentum direction
- WASM-compatible math kernel (`cargo build --target wasm32-unknown-unknown --lib`)

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
use lambert_izzo::{lambert, TransferWay};
use nalgebra::Vector3;

// LEO → MEO Hohmann transfer.
let mu_km3_s2 = 398_600.441_8;
let r1_km = Vector3::new(7000.0, 0.0, 0.0);
let r2_km = Vector3::new(-12_000.0, 1.0, 0.0); // 1 km off-axis avoids colinearity
let a_km = (7000.0 + 12_000.0) / 2.0;
let tof_s = std::f64::consts::PI * (a_km.powi(3) / mu_km3_s2).sqrt();

let solutions = lambert(r1_km, r2_km, tof_s, mu_km3_s2, TransferWay::Short, 0).unwrap();
let sol = &solutions[0];
println!("v1 = {} km/s", sol.v1_km_s);
println!("v2 = {} km/s", sol.v2_km_s);
println!("converged in {} iterations", sol.iters);
```

The signature:

```rust
pub fn lambert(
    r1_km: Vector3<f64>,    // initial position (km), any inertial frame
    r2_km: Vector3<f64>,    // final position (km), same frame
    tof_s: f64,             // time of flight (s), > 0
    mu_km3_s2: f64,         // gravitational parameter (km³/s²), > 0
    way: TransferWay,       // Short or Long way around the transfer plane
    multi_revs: u32,        // max revolution count (0 = single-rev only)
) -> Result<Vec<LambertSolution>, LambertError>
```

For `multi_revs = N`, you get up to `1 + 2 · N_max` solutions, where
`N_max = min(multi_revs, ⌊T/π⌋)`. Single-rev is always index 0; multi-rev
branches (if any) follow as `(long-period, short-period)` pairs.

`LambertError` is a `thiserror` enum with structured fields:

```rust
match lambert(r1_km, r2_km, tof_s, mu_km3_s2, TransferWay::Short, 0) {
    Ok(sols) => /* … */,
    Err(LambertError::NonPositiveTimeOfFlight { tof_s })       => /* … */,
    Err(LambertError::NonPositiveMu { mu_km3_s2 })             => /* … */,
    Err(LambertError::DegeneratePositionVector { which, norm_km }) => /* … */,
    Err(LambertError::CollinearGeometry { sin_angle })         => /* … */,
    Err(LambertError::NoConvergence { iterations, last_step, n_revs }) => /* … */,
    Err(_) => /* … */,
}
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
| Single-rev | 100%        | 2.083     | 2.1       | 3         | 1.18e-11 | 3.16e-12 |
| Multi-rev  | 100%        | 2.992     | 3.3       | 6         | 2.77e-14 | 1.14e-13 |

Iteration counts match the paper's reported figures.

## Building

```
cargo build --release
cargo test --release
cargo run --release --example demo
cargo run --release --example stress
```

Toolchain pinned via `rust-toolchain.toml` (1.88.0). Edition 2024.
Dependencies: [`nalgebra`](https://nalgebra.org) for `Vector3<f64>` arithmetic
and [`thiserror`](https://docs.rs/thiserror) for the error type.

## Implementation notes

- Modular layout under `src/`:
  - `constants.rs` — every named tolerance with rationale.
  - `error.rs` — structured `LambertError` enum.
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

## License

MIT OR Apache-2.0
