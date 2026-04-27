# Lambert's problem — a fast intro

This crate solves Lambert's two-point boundary value problem under
two-body gravity. If you've never encountered the problem, this page is
the on-ramp before the [paper](izzo.pdf).

## What is Lambert's problem?

Given:

- A starting position `r1` (in some inertial frame).
- A target position `r2` (in the same frame).
- A time of flight `tof` between the two.
- A gravitational parameter `μ` of the central body (Earth, Sun, …).

Find: the velocity `v1` at `r1` such that the resulting Keplerian orbit
arrives at `r2` exactly `tof` seconds later. Plus the velocity `v2` at
arrival.

That's the whole problem. It's the workhorse of mission design — every
porkchop plot, rendezvous targeting, and intercept calculation reduces to
"pick `(departure_time, arrival_time)` and solve a Lambert problem."

## Why is it non-trivial?

Two reasons.

**One.** The relationship between `tof` and the orbital geometry is
transcendental — there's no closed-form `v1 = f(r1, r2, tof, μ)`. The
solver has to root-find on an auxiliary scalar (Izzo calls it `x`,
following Lancaster and Blanchard) until the numerically-evaluated
time-of-flight matches the requested `tof`.

**Two.** The same boundary problem can have *multiple* valid orbits.
Specifically:

- Two **transfer directions**: short-way (the geodesic, `θ ≤ π`) and
  long-way (`θ > π`). Same start, same end — different ways around.
- Multiple **revolution counts**: for long enough `tof`, the spacecraft
  can wrap around the central body once, twice, … `M` times before
  arriving. Each revolution count `M ≥ 1` admits *two* trajectories: a
  long-period one (more time near apoapsis) and a short-period one
  (more time near periapsis).

This crate returns all of them, as a typed `LambertSolutions { single,
multi: [(M=1 long, M=1 short), (M=2 long, M=2 short), ...] }`.

## Picturing the geometry

ASCII sketch — Earth at the origin, `r1` on the +X axis, `r2` 90° out on
the +Y axis:

```
                     r2
                     |
                     ●
                    /|
                   / |
                  /  |
                 /   |
       short ●— /    |
       way  ↗  /     |
              /      |
             /       |
            ●━━━━━━━━━●━━━━━━━ r1
            ↘  long way (the rest of the way around)
```

Both arcs end at `r2`, but one is ~270° while the other is ~90°. That
choice is the `TransferWay::Short` / `TransferWay::Long` argument.

The angular-momentum direction (whether the orbit goes counter-clockwise
or clockwise as seen from above) is set by the *order* of `r1` and
`r2` — `r1 × r2` defines the orbit plane's normal. Swap them to flip
prograde/retrograde.

## The Lancaster–Blanchard `x`

Internally, the solver works in a single dimensionless parameter `x ∈
(-∞, ∞)`. Roughly:

- `x = -1`: parabolic escape from the long side.
- `x = 0`: a "balanced" elliptic transfer.
- `x = 1`: parabolic escape from the short side.
- `|x| > 1`: hyperbolic.

For each `(r1, r2, μ, way)`, the time of flight `T(x)` is a smooth-ish
function with a minimum and a single root for any feasible `tof`.
Multi-rev branches add `M·π` periods, splitting `T(x)` into pieces with
their own minima — hence the long-period and short-period roots per `M`.

Returned in `SolverDiagnostics::lancaster_blanchard_x` if you want it.
The long-period branch always has the smaller `x` for a given `M`.

## What the three TOF regimes are about

Numerically, the time-of-flight formula for a given `x` is unstable near
`x = 1` (the parabolic point). The solver dispatches:

- **Battin's hypergeometric series** (Izzo Eq. 20) for `|x − 1| ≤ 0.01`
  — exact near the singularity.
- **Lancaster–Blanchard form** (Izzo Eq. 18) for `0.01 < |x − 1| ≤ 0.2`
  — clean middle band.
- **Lagrange form** (Izzo Eq. 9) for `|x − 1| > 0.2` — stable far from
  parabolic.

You don't pick the regime — `tof::x_to_tof` does. The thresholds live
in `constants.rs` with rationale.

## What this crate doesn't do

- Escape trajectories from one body to another (e.g. Earth-to-Mars
  *patched-conic* with sphere-of-influence handoffs). Lambert solves
  the inertial transfer; you assemble the patched-conic story around it.
- Lunar swing-by or higher-order n-body effects. Pure two-body.
- Low-thrust transfers. Lambert is impulsive — instant Δv at endpoints.
- Constraint-satisfaction (e.g. "land between these two delta-V budgets"
  — that's an outer optimization loop calling Lambert in the inner
  loop).

## Further reading

- Izzo, D. (2014). [*Revisiting Lambert's
  problem*](https://arxiv.org/abs/1403.2705). The reference algorithm
  this crate implements; included as `docs/izzo.pdf`.
- Battin, R. H. (1999). *An Introduction to the Mathematics and Methods
  of Astrodynamics*, Revised Edition. Chapter 7 covers the transfer
  problem and the hypergeometric series form used here.
- Lancaster, E. R. and Blanchard, R. C. (1969). *A Unified Form of
  Lambert's Theorem*. NASA TN D-5368. The origin of the dimensionless
  `x` parameter and the canonical form of `T(x)`.
