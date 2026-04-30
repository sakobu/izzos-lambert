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

A call to `lambert(&input)` runs three top-level stages — input
validation + geometry pre-compute (`from_inputs`), the `find_xy` core
(which itself nests multi-rev enumeration → Householder iteration → TOF
evaluation), and velocity reconstruction (`build_solutions`).

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
┌──── find_xy: branch loop ────────────────────────────────────┐
│  3. Multi-rev enumeration  → root_finding.rs (find_xy loop)  │
│  4. Root finding (per M)   → root_finding.rs (Householder +  │
│                              Halley T_min)                   │
│  5. TOF evaluation (in 4)  → tof.rs (Battin / Lancaster /    │
│                              Lagrange)                       │
└──────────────────────────────────────────────────────────────┘
     │
     ▼
6. Velocity reconstruction → lib.rs (reconstruct → build_solutions)
     │
     ▼
LambertSolutions
```

The driver is `lambert` in
[`crates/lambert_izzo/src/lib.rs`](crates/lambert_izzo/src/lib.rs).
Three lines:

```rust
let geom = Geometry::from_inputs(input.r1, input.r2, input.tof, input.mu, input.way)?;
let roots = find_xy(&geom, input.revolutions)?;
Ok(build_solutions(&geom, &roots))
```

## 1. Input validation

Public errors are defined in `error.rs`; the validation logic itself
runs at the top of `Geometry::from_inputs`. Each check returns a typed
`LambertError` variant with structured field data — pattern-match on
the value, never parse a string.

- `LambertError` enum:
  [`crates/lambert_izzo/src/error.rs`](crates/lambert_izzo/src/error.rs).
  Variants: `NonFiniteInput`, `NonPositiveTimeOfFlight`,
  `NonPositiveMu`, `DegeneratePositionVector`, `CollinearGeometry`,
  `NoConvergence`, `SingularDenominator`.
- Endpoint and parameter identity (so callers can pattern-match
  rather than read messages): `Position::{R1, R2}` and
  `NonFiniteParameter`, both in
  [`error.rs`](crates/lambert_izzo/src/error.rs).
- Construction-time bound on revolution counts (a separate error type,
  raised before the solver is even called): `RevsOutOfRange` in
  [`bounded_revs.rs`](crates/lambert_izzo/src/bounded_revs.rs).
- Validation site: the top of `Geometry::from_inputs` in
  [`geometry.rs`](crates/lambert_izzo/src/geometry.rs) — finiteness,
  sign, and norm-floor checks; collinearity is checked once the cross
  product is available.

## 2. Geometry pre-computation

`Geometry::from_inputs` is the only constructor; the resulting
`Geometry` struct is then threaded read-only through stages 3–6 so the
chord, semi-perimeter, and λ are never re-derived in two places.

- Struct definition: `Geometry` in
  [`geometry.rs`](crates/lambert_izzo/src/geometry.rs). Fields
  are paper symbols (`lambda`, `big_t`, `gamma`, `rho`, `sigma`, plus
  the unit-vector basis `ir1`, `ir2`, `it1`, `it2` and norms `r1n`,
  `r2n`).
- Build site: `Geometry::from_inputs` in
  [`geometry.rs`](crates/lambert_izzo/src/geometry.rs).
- λ sign convention (Izzo Eq. 7) and short/long-way handling: inside
  `from_inputs`.
- Round-off guards — both apply `.max(0.0)` before a `sqrt` to absorb
  the few-ulp negatives that appear when `c` rounds slightly above
  `s`, or `ρ` slightly above `1`.

## 3. Multi-rev branch enumeration

The Izzo formulation admits the always-present single-revolution
solution plus up to `⌊T/π⌋` multi-revolution branches, each of which
yields a long-period and a short-period root. The crate caps
enumeration at `BoundedRevs::MAX = 32`, type-enforced at construction
time so the bounded return collection always fits without truncation.

- `RevolutionBudget` enum (`SingleOnly` / `UpTo(BoundedRevs)`):
  [`lib.rs`](crates/lambert_izzo/src/lib.rs).
- `BoundedRevs` newtype + `MAX = 32`:
  [`bounded_revs.rs`](crates/lambert_izzo/src/bounded_revs.rs).
- Iterator over `1..=M` as validated `BoundedRevs` values, used to
  drive the loop in stage 4: `RevolutionBudget::iter_revs` in
  [`lib.rs`](crates/lambert_izzo/src/lib.rs).
- Branch loop site: inside `find_xy` in
  [`root_finding.rs`](crates/lambert_izzo/src/root_finding.rs).
- **Silent-skip semantics** (the load-bearing detail): higher-`M`
  branches are dropped from the returned `multi` set when
  `T_min(M) > tof`, and the loop stops at the first infeasible `M`
  (since `T_min` is monotone in `M`). Two checks gate this:
  - Quick reject `tof < M·π`: every branch needs at least `M·π`
    non-dimensional time.
  - Boundary-zone Halley `T_min` confirmation when `tof` lies below
    the analytic minimum at `x = 0`.
- Caller-side observability — the highest `M` actually solved:
  `LambertSolutions::max_feasible_revs` in
  [`lib.rs`](crates/lambert_izzo/src/lib.rs).

## 4. Root finding

Householder's third-order method on `T(x) − T_target = 0`, started
from analytic guesses (Izzo Eq. 30 single-rev, Eq. 31 multi-rev). The
iteration consumes the analytic derivatives of `T(x)` (Eq. 22)
produced by stage 5 — see `tof_derivatives_with_y` below.

> **Eq. 30 typo caveat.** The paper's Eq. 30 displays the middle-branch
> initial guess as `(T0/T)^(log2(T1/T0)) − 1`, which fails the boundary
> `x = 1` at `T = T1`. The textual derivation just before Eq. 30
> (page 12) gives the correct form `(T0/T)^(ln 2 / ln(T0/T1)) − 1`, and
> the reference PyKEP implementation uses the same corrected form.
> `lambert_izzo` follows that corrected form — see the in-source
> comment in `initial_guess_single_rev`.

- Driver: `find_xy` in
  [`root_finding.rs`](crates/lambert_izzo/src/root_finding.rs).
  Returns the always-present single-rev `Root` plus an
  `ArrayVec<RootPair, MAX_MULTI_REV_PAIRS>` of feasible multi-rev pairs.
- Single-rev initial guess (Eq. 30 plus the derivative-matched
  hyperbolic starter for the `T ≤ T1` branch):
  `initial_guess_single_rev` in
  [`root_finding.rs`](crates/lambert_izzo/src/root_finding.rs).
- Multi-rev initial guesses (Eq. 31, returns the long-period and
  short-period asymptote starters as `(x0l, x0r)`):
  `initial_guess_multi_rev` in
  [`root_finding.rs`](crates/lambert_izzo/src/root_finding.rs).
- Householder iteration (Eq. 22):
  [`root_finding.rs`](crates/lambert_izzo/src/root_finding.rs).
  Returns `LambertError::SingularDenominator` if the denominator
  vanishes algebraically, or `LambertError::NoConvergence` if the
  `HOUSEHOLDER_MAX_ITERS = 15` cap is reached.
- Halley iteration on `dT/dx = 0` for the multi-rev `T_min` feasibility
  check used by stage 3: `halley_t_min` in
  [`root_finding.rs`](crates/lambert_izzo/src/root_finding.rs).
  Cubic convergence reaches `f64` precision well inside
  `HALLEY_MAX_ITERS = 12`.
- All convergence tolerances and iteration caps:
  [`constants.rs`](crates/lambert_izzo/src/constants.rs).

## 5. TOF evaluation

`T(x, λ, M)` is evaluated in three regimes that blend at fixed
thresholds in `|x − 1|`. The dispatcher picks one; the Householder
loop never sees the choice. For the *why* behind the regimes (which
formula loses precision where, and what the thresholds are tuned
against), read the matching section in
[`docs/concepts.md`](docs/concepts.md).

- Dispatcher: `x_to_tof_with_y` in
  [`crates/lambert_izzo/src/tof.rs`](crates/lambert_izzo/src/tof.rs).
- Battin hypergeometric series (Eq. 20), `|x − 1| ≤ BATTIN_THRESHOLD`,
  single-rev only.
- Lancaster–Blanchard (Eq. 18),
  `BATTIN_THRESHOLD < |x − 1| ≤ LAGRANGE_THRESHOLD`.
- Lagrange (Eq. 9), elsewhere.
- Threshold constants: `BATTIN_THRESHOLD`, `LAGRANGE_THRESHOLD` in
  [`constants.rs`](crates/lambert_izzo/src/constants.rs).
- Analytic derivatives `(dT/dx, d²T/dx², d³T/dx³)` (Eq. 22) consumed
  by Householder in stage 4: `tof_derivatives_with_y` in
  [`tof.rs`](crates/lambert_izzo/src/tof.rs).

## 6. Velocity reconstruction

Each converged `(x, y)` root becomes a `(v1, v2)` pair via Izzo
Algorithm 1 — radial and tangential components computed from λ, ρ, σ,
γ, then projected onto the in-plane basis built in stage 2.

- Per-root reconstruction: `reconstruct` in
  [`crates/lambert_izzo/src/lib.rs`](crates/lambert_izzo/src/lib.rs).
- Driver wrapping `reconstruct` over the single-rev root and every
  multi-rev pair, plus the parallel `LambertDiagnostics` assembly:
  `build_solutions` in
  [`lib.rs`](crates/lambert_izzo/src/lib.rs).
- The output type and its diagnostics counterpart: `LambertSolutions`
  (with `max_feasible_revs`) and `LambertDiagnostics`, both in
  [`lib.rs`](crates/lambert_izzo/src/lib.rs).

## See also

- [`docs/concepts.md`](docs/concepts.md) — what Lambert's problem is,
  the role of the dimensionless `x`, and why three TOF regimes exist.
- [`docs/izzo.pdf`](docs/izzo.pdf) — D. Izzo (2014), the algorithm
  reference. Inline comments cite this paper as `Eq. N` / `Algorithm N`.
