# Architecture

A code map of `lambert_izzo`: how a single call to `lambert(&input)`
flows through the crate, and where each stage lives in the source.

The crate solves Lambert's two-point boundary-value problem under
two-body gravity using D. Izzo's revisited algorithm
([arXiv:1403.2705](https://arxiv.org/abs/1403.2705); PDF at
[`docs/izzo.pdf`](docs/izzo.pdf)). For an intro to *what* the problem
is and *why* three TOF regimes exist, read
[`docs/concepts.md`](docs/concepts.md) first — this document picks up
where that one leaves off and points at code.

The public surface is plain `f64` and `[f64; 3]` — no math-library
dependency. The math is dimensionally homogeneous and frame-invariant:
pass `r1` / `r2` in any consistent inertial frame and consistent units,
and the returned velocities come back in that same frame and units.

## Pipeline overview

A call to `lambert(&input)` runs six stages in order. Each has its own
section below.

```
LambertInput
     │
     ▼
1. Input validation       → geometry.rs (Geometry::from_inputs)
     │
     ▼
2. Geometry pre-compute   → geometry.rs (lambda, big_t, gamma, ρ, σ, …)
     │
     ▼
3. Multi-rev enumeration  → root_finding.rs (find_xy, branch loop)
     │
     ▼
4. Root finding           → root_finding.rs (Householder + Halley T_min)
     │
     ▼
5. TOF evaluation         → tof.rs (Battin / Lancaster / Lagrange)
     │
     ▼
6. Velocity reconstruction → lib.rs (reconstruct → build_solutions)
     │
     ▼
LambertSolutions
```

The driver is `lambert` at
[`crates/lambert_izzo/src/lib.rs:454`](crates/lambert_izzo/src/lib.rs).
Three lines:

```rust
let geom = Geometry::from_inputs(input.r1, input.r2, input.tof, input.mu, input.way)?;
let roots = find_xy(&geom, input.revolutions)?;
Ok(build_solutions(&geom, &roots))
```

Stage 1 + 2 happen inside `from_inputs`, stages 3 + 4 + 5 inside
`find_xy`, and stage 6 inside `build_solutions`.

## 1. Input validation

Public errors are defined in `error.rs`; the validation logic itself
runs at the top of `Geometry::from_inputs`. Each check returns a typed
`LambertError` variant with structured field data — pattern-match on
the value, never parse a string.

- `LambertError` enum:
  [`crates/lambert_izzo/src/error.rs:80-149`](crates/lambert_izzo/src/error.rs).
  Variants: `NonFiniteInput`, `NonPositiveTimeOfFlight`,
  `NonPositiveMu`, `DegeneratePositionVector`, `CollinearGeometry`,
  `NoConvergence`, `SingularDenominator`.
- Endpoint and parameter identity (so callers can pattern-match
  rather than read messages):
  [`error.rs:6-20`](crates/lambert_izzo/src/error.rs) (`Position::{R1, R2}`),
  [`error.rs:29-46`](crates/lambert_izzo/src/error.rs) (`NonFiniteParameter`).
- Construction-time bound on revolution counts (a separate error type,
  raised before the solver is even called):
  [`crates/lambert_izzo/src/bounded_revs.rs:91-99`](crates/lambert_izzo/src/bounded_revs.rs)
  (`RevsOutOfRange`).
- Validation site:
  [`geometry.rs:64-96`](crates/lambert_izzo/src/geometry.rs) — finiteness,
  sign, and norm-floor checks; collinearity is checked at line 106 once
  the cross product is available.

## 2. Geometry pre-computation

`Geometry::from_inputs` is the only constructor; the resulting
`Geometry` struct is then threaded read-only through stages 3–6 so the
chord, semi-perimeter, and λ are never re-derived in two places.

- Struct definition:
  [`geometry.rs:24-47`](crates/lambert_izzo/src/geometry.rs). Fields
  are paper symbols (`lambda`, `big_t`, `gamma`, `rho`, `sigma`, plus
  the unit-vector basis `ir1`, `ir2`, `it1`, `it2` and norms `r1n`,
  `r2n`).
- Build site:
  [`geometry.rs:64-144`](crates/lambert_izzo/src/geometry.rs)
  (`Geometry::from_inputs`).
- λ sign convention (Izzo Eq. 7) and short/long-way handling:
  [`geometry.rs:111-120`](crates/lambert_izzo/src/geometry.rs).
- Round-off guards — both apply `.max(0.0)` before a `sqrt` to absorb
  the few-ulp negatives that appear when `c` rounds slightly above
  `s`, or `ρ` slightly above `1`:
  [`geometry.rs:114`](crates/lambert_izzo/src/geometry.rs) (λ),
  [`geometry.rs:129`](crates/lambert_izzo/src/geometry.rs) (σ).

## 3. Multi-rev branch enumeration

The Izzo formulation admits the always-present single-revolution
solution plus up to `⌊T/π⌋` multi-revolution branches, each of which
yields a long-period and a short-period root. The crate caps
enumeration at `BoundedRevs::MAX = 32`, type-enforced at construction
time so the bounded return collection always fits without truncation.

- `RevolutionBudget` enum (`SingleOnly` / `UpTo(BoundedRevs)`):
  [`lib.rs:146-212`](crates/lambert_izzo/src/lib.rs).
- `BoundedRevs` newtype + `MAX = 32`:
  [`bounded_revs.rs:14-68`](crates/lambert_izzo/src/bounded_revs.rs).
- Iterator over `1..=M` as validated `BoundedRevs` values, used to
  drive the loop in stage 4:
  [`lib.rs:207-211`](crates/lambert_izzo/src/lib.rs)
  (`RevolutionBudget::iter_revs`).
- Branch loop site:
  [`crates/lambert_izzo/src/root_finding.rs:95-137`](crates/lambert_izzo/src/root_finding.rs)
  inside `find_xy`.
- **Silent-skip semantics** (the load-bearing detail): higher-`M`
  branches are dropped from the returned `multi` set when
  `T_min(M) > tof`, and the loop stops at the first infeasible `M`
  (since `T_min` is monotone in `M`). Two checks gate this:
  - Quick reject `tof < M·π`: every branch needs at least `M·π`
    non-dimensional time
    ([`root_finding.rs:99-102`](crates/lambert_izzo/src/root_finding.rs)).
  - Boundary-zone Halley `T_min` confirmation when `tof` lies below
    the analytic minimum at `x = 0`
    ([`root_finding.rs:104-112`](crates/lambert_izzo/src/root_finding.rs)).
- Caller-side observability — the highest `M` actually solved:
  [`lib.rs:297-299`](crates/lambert_izzo/src/lib.rs)
  (`LambertSolutions::max_feasible_revs`).

## 4. Root finding

Householder's third-order method on `T(x) − T_target = 0`, started
from analytic guesses (Izzo Eq. 30 single-rev, Eq. 31 multi-rev). The
iteration consumes the analytic derivatives of `T(x)` (Eq. 22)
produced by stage 5 — see `tof_derivatives_with_y` below.

- Driver: [`root_finding.rs:65-139`](crates/lambert_izzo/src/root_finding.rs)
  (`find_xy`). Returns the always-present single-rev `Root` plus an
  `ArrayVec<RootPair, MAX_MULTI_REV_PAIRS>` of feasible multi-rev pairs.
- Single-rev initial guess (Eq. 30 plus the derivative-matched
  hyperbolic starter for the `T ≤ T1` branch):
  [`root_finding.rs:143-158`](crates/lambert_izzo/src/root_finding.rs).
- Multi-rev initial guesses (Eq. 31, returns the long-period and
  short-period asymptote starters as `(x0l, x0r)`):
  [`root_finding.rs:163-170`](crates/lambert_izzo/src/root_finding.rs).
- Householder iteration (Eq. 22):
  [`root_finding.rs:174-202`](crates/lambert_izzo/src/root_finding.rs).
  Returns `LambertError::SingularDenominator` if the denominator
  vanishes algebraically, or `LambertError::NoConvergence` if the
  `HOUSEHOLDER_MAX_ITERS = 15` cap is reached.
- Halley iteration on `dT/dx = 0` for the multi-rev `T_min` feasibility
  check used by stage 3:
  [`root_finding.rs:215-240`](crates/lambert_izzo/src/root_finding.rs)
  (`halley_t_min`). Cubic convergence reaches `f64` precision well
  inside `HALLEY_MAX_ITERS = 12`.
- All convergence tolerances and iteration caps:
  [`crates/lambert_izzo/src/constants.rs:25-55`](crates/lambert_izzo/src/constants.rs).

## 5. TOF evaluation

`T(x, λ, M)` is evaluated in three regimes that blend at fixed
thresholds in `|x − 1|`. The dispatcher picks one; the Householder
loop never sees the choice. For the *why* behind the regimes (which
formula loses precision where, and what the thresholds are tuned
against), read the matching section in
[`docs/concepts.md`](docs/concepts.md).

- Dispatcher:
  [`crates/lambert_izzo/src/tof.rs:26-35`](crates/lambert_izzo/src/tof.rs)
  (`x_to_tof_with_y`).
- Battin hypergeometric series (Eq. 20), `|x − 1| ≤ 0.01`,
  single-rev only:
  [`tof.rs:96-101`](crates/lambert_izzo/src/tof.rs).
- Lancaster–Blanchard (Eq. 18), `0.01 < |x − 1| ≤ 0.2`:
  [`tof.rs:73-77`](crates/lambert_izzo/src/tof.rs).
- Lagrange (Eq. 9), elsewhere:
  [`tof.rs:49-69`](crates/lambert_izzo/src/tof.rs).
- Threshold values:
  [`constants.rs:68`](crates/lambert_izzo/src/constants.rs)
  (`BATTIN_THRESHOLD = 0.01`),
  [`constants.rs:74`](crates/lambert_izzo/src/constants.rs)
  (`LAGRANGE_THRESHOLD = 0.2`).
- Analytic derivatives `(dT/dx, d²T/dx², d³T/dx³)` (Eq. 22) consumed
  by Householder in stage 4:
  [`tof.rs:127-142`](crates/lambert_izzo/src/tof.rs)
  (`tof_derivatives_with_y`).

## 6. Velocity reconstruction

Each converged `(x, y)` root becomes a `(v1, v2)` pair via Izzo
Algorithm 1 — radial and tangential components computed from λ, ρ, σ,
γ, then projected onto the in-plane basis built in stage 2.

- Per-root reconstruction:
  [`crates/lambert_izzo/src/lib.rs:489-514`](crates/lambert_izzo/src/lib.rs)
  (`reconstruct`).
- Driver wrapping `reconstruct` over the single-rev root and every
  multi-rev pair, plus the parallel `LambertDiagnostics` assembly:
  [`lib.rs:516-550`](crates/lambert_izzo/src/lib.rs)
  (`build_solutions`).
- The output type and its diagnostics counterpart:
  [`lib.rs:271-300`](crates/lambert_izzo/src/lib.rs)
  (`LambertSolutions`, `max_feasible_revs`),
  [`lib.rs:326-334`](crates/lambert_izzo/src/lib.rs)
  (`LambertDiagnostics`).

## See also

- [`docs/concepts.md`](docs/concepts.md) — what Lambert's problem is,
  the role of the dimensionless `x`, and why three TOF regimes exist.
- [`docs/izzo.pdf`](docs/izzo.pdf) — D. Izzo (2014), the algorithm
  reference. Inline comments cite this paper as `Eq. N` / `Algorithm N`.
