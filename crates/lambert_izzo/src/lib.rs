//! Izzo's revisited Lambert solver — single + multi-revolution, short/long way.
//!
//! Reference: D. Izzo, *Revisiting Lambert's problem*, Celestial Mechanics &
//! Dynamical Astronomy, 2014. arXiv:1403.2705. PDF in `docs/izzo.pdf`.
//!
//! Inline `Eq. N` / `Algorithm N` references in the source point to that paper.
//!
//! # Units
//!
//! The crate is unit-agnostic at the type level — all quantities are plain
//! `f64` and `[f64; 3]`. The convention used throughout the docs and
//! examples is SI for astrodynamics work:
//!
//! | Quantity                | Unit   |
//! | ----------------------- | ------ |
//! | Position                | km     |
//! | Velocity                | km/s   |
//! | Time of flight          | s      |
//! | Gravitational parameter | km³/s² |
//!
//! Any consistent unit system works — pass `r` in meters and `mu` in m³/s²
//! and you get velocities in m/s. The math is dimensionally homogeneous.
//!
//! The algorithm is also frame-invariant under any inertial frame — pass
//! `r1` and `r2` in the same inertial frame (ECI for Earth orbits, HCRS
//! for solar transfers, etc.) and the returned velocities are in that
//! same frame. The function signature is frame-agnostic; the calling
//! code's variable names carry the frame info.
//!
//! Position and velocity vectors are plain `[f64; 3]` arrays; the crate has
//! no hard math-library dependency. Both [`nalgebra::Vector3<f64>`] and
//! [`glam::DVec3`] already convert to/from `[f64; 3]` natively, so callers
//! using either library can pass and receive vectors without an explicit
//! interop layer.
//!
//! ```ignore
//! // nalgebra:
//! let r1: [f64; 3] = nalgebra::Vector3::new(7000.0, 0.0, 0.0).into();
//! let v1: nalgebra::Vector3<f64> = solution.single.v1.into();
//!
//! // glam:
//! let r2 = glam::DVec3::new(0.0, 7000.0, 0.0).to_array();
//! let v2 = glam::DVec3::from_array(solution.single.v2);
//! ```
//!
//! # Example
//!
//! ```
//! use lambert_izzo::{lambert, LambertError, RevolutionBudget, TransferWay};
//!
//! # fn main() -> Result<(), LambertError> {
//! // LEO → LEO 90° transfer at 7000 km altitude.
//! let mu = 398_600.4418;
//! let r1 = [7000.0, 0.0, 0.0];
//! let r2 = [0.0, 7000.0, 0.0];
//! let tof = core::f64::consts::PI / 2.0 * (7000.0_f64.powi(3) / mu).sqrt();
//!
//! let solutions = lambert(
//!     r1, r2, tof, mu,
//!     TransferWay::Short, RevolutionBudget::SingleOnly,
//! )?;
//! assert!(solutions.multi.is_empty());
//! let v1 = solutions.single.v1;
//! # let _ = v1;
//! # Ok(())
//! # }
//! ```

#![cfg_attr(not(test), no_std)]
#![warn(missing_docs)]
#![warn(clippy::pedantic)]
#![warn(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::unreachable
)]
#![allow(clippy::module_name_repetitions)] // LambertError, LambertSolution, LambertSolutions, LambertDiagnostics

use arrayvec::ArrayVec;

/// Hard upper bound on the number of multi-revolution pairs returned.
///
/// The Izzo formulation admits up to `⌊T/π⌋` multi-rev branches, which can
/// be arbitrarily large for very long times of flight — but practical
/// missions almost never exceed `M = 5`. This cap keeps the bounded return
/// type a fixed stack size.
///
/// Callers passing [`RevolutionBudget::up_to`] above this cap get the
/// truncated set rather than an error.
pub const MAX_MULTI_REV_PAIRS: usize = 32;

/// Direction around the transfer plane from `r1` to `r2`.
///
/// `Short` is the geodesic transfer (`θ ≤ π`); `Long` traverses the other
/// way (`θ > π`). This is independent of prograde/retrograde — the orbit's
/// angular-momentum direction is set by the order of the `(r1, r2)` arguments.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum TransferWay {
    /// `θ ≤ π` — the short geodesic arc.
    Short,
    /// `θ > π` — the long way around the transfer plane.
    Long,
}

/// Maximum number of complete revolutions to consider beyond single-rev.
///
/// Multi-revolution branches admit two solutions per revolution count `M`
/// (long-period and short-period), so the total solution count is
/// `1 + 2 · min(max(), ⌊T/π⌋)` adjusted downward when a branch's `T_min`
/// exceeds the requested time of flight, and clamped at
/// [`MAX_MULTI_REV_PAIRS`] regardless of the requested value.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum RevolutionBudget {
    /// Solve the single-revolution case only — always exactly one solution.
    #[default]
    SingleOnly,
    /// Search up to `M` complete revolutions inclusive (`M ≥ 1`).
    UpTo(core::num::NonZeroU32),
}

impl RevolutionBudget {
    /// Convenience constructor: `up_to(0)` collapses to `SingleOnly`.
    #[must_use]
    pub fn up_to(m: u32) -> Self {
        match core::num::NonZeroU32::new(m) {
            Some(nz) => Self::UpTo(nz),
            None => Self::SingleOnly,
        }
    }

    /// The maximum revolution count this budget will search (`0` for single-only).
    #[must_use]
    pub fn max(self) -> u32 {
        match self {
            Self::SingleOnly => 0,
            Self::UpTo(n) => n.get(),
        }
    }
}

pub mod constants;
mod error;
mod geometry;
mod root_finding;
mod tof;
mod vec3;

pub use error::{LambertError, NonFiniteParameter, Position};

use geometry::Geometry;
use root_finding::{Root, find_xy};

/// One Lambert transfer trajectory.
///
/// Pure trajectory data — start and end velocities. Solver diagnostics
/// (iteration count, Lancaster–Blanchard `x`) are kept out of this type so
/// the common API path is lean; use [`solve_with_diagnostics`] when you
/// need them.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct LambertSolution {
    /// Velocity at `r1` (same units and frame as the inputs).
    pub v1: [f64; 3],
    /// Velocity at `r2` (same units and frame as the inputs).
    pub v2: [f64; 3],
}

/// One multi-revolution pair: long-period and short-period trajectories
/// for a given revolution count.
///
/// The Izzo formulation admits exactly two trajectories for each `M ≥ 1`:
/// the long-period branch (smaller `x`, more time near apoapsis) and the
/// short-period branch (larger `x`, more time near periapsis).
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct MultiRevPair {
    /// Branch revolution count (`>= 1`).
    pub n_revs: u32,
    /// Long-period trajectory (smaller Lancaster–Blanchard `x`).
    pub long_period: LambertSolution,
    /// Short-period trajectory (larger Lancaster–Blanchard `x`).
    pub short_period: LambertSolution,
}

/// All Lambert trajectories for a given boundary problem and revolution budget.
///
/// Always carries the single-revolution trajectory; `multi` lists every
/// reachable multi-rev pair in ascending `M` order, capped at
/// [`MAX_MULTI_REV_PAIRS`].
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct LambertSolutions {
    /// Single-revolution trajectory — always present.
    pub single: LambertSolution,
    /// Multi-revolution pairs in ascending `M` order; empty for
    /// [`RevolutionBudget::SingleOnly`] or when no multi-rev branches are
    /// feasible for the given time of flight.
    pub multi: ArrayVec<MultiRevPair, MAX_MULTI_REV_PAIRS>,
}

/// Solutions for both transfer ways, computed from the same boundary inputs.
///
/// Returned by [`lambert_both_ways`]. The two halves are independent —
/// either one may be empty in `multi` while the other is populated, since
/// `T_min(M)` differs between the short and long forms.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct BothWaysSolutions {
    /// Short-way trajectories (`θ ≤ π`).
    pub short: LambertSolutions,
    /// Long-way trajectories (`θ > π`).
    pub long: LambertSolutions,
}

/// Diagnostic data for one converged Householder solve.
///
/// Useful for debugging or distinguishing multi-rev branches; not part of
/// the trajectory answer. Returned by [`solve_with_diagnostics`].
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct SolverDiagnostics {
    /// Householder iterations used to converge.
    pub iters: u32,
    /// Final value of Izzo's free parameter (Lancaster–Blanchard variable),
    /// dimensionless. For a given `M`, the long-period branch has the
    /// smaller value, the short-period branch the larger.
    pub lancaster_blanchard_x: f64,
}

/// Diagnostics for one multi-rev pair.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct MultiRevPairDiagnostics {
    /// Branch revolution count (`>= 1`).
    pub n_revs: u32,
    /// Long-period branch diagnostics.
    pub long_period: SolverDiagnostics,
    /// Short-period branch diagnostics.
    pub short_period: SolverDiagnostics,
}

/// Diagnostics structure mirroring [`LambertSolutions`].
///
/// Returned alongside the solutions by [`solve_with_diagnostics`].
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct LambertDiagnostics {
    /// Single-rev solver diagnostics.
    pub single: SolverDiagnostics,
    /// Multi-rev pair diagnostics in the same order as `LambertSolutions::multi`.
    pub multi: ArrayVec<MultiRevPairDiagnostics, MAX_MULTI_REV_PAIRS>,
}

/// One Lambert call's inputs, packaged for batch processing.
///
/// Same fields as [`lambert`]'s parameter list, in struct form so callers
/// can build a slice (e.g. for porkchop plots) and stream it through
/// [`lambert_iter`] or [`lambert_par_iter`].
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct LambertInput {
    /// Initial position, any consistent inertial frame.
    pub r1: [f64; 3],
    /// Final position, same frame as `r1`.
    pub r2: [f64; 3],
    /// Time of flight, `> 0`.
    pub tof: f64,
    /// Gravitational parameter of the central body, `> 0`.
    pub mu: f64,
    /// Short or long way around the transfer plane.
    pub way: TransferWay,
    /// Revolution budget — see [`RevolutionBudget`].
    pub revolutions: RevolutionBudget,
}

impl LambertInput {
    /// Solve this single input — convenience wrapper around [`lambert`].
    ///
    /// # Errors
    ///
    /// Same conditions as [`lambert`].
    pub fn solve(self) -> Result<LambertSolutions, LambertError> {
        lambert(self.r1, self.r2, self.tof, self.mu, self.way, self.revolutions)
    }
}

/// Sequential batch iterator over Lambert inputs.
///
/// Allocation-free; just maps over the input slice. Useful for
/// porkchop-plot-style workloads where the caller computes one Lambert
/// solution per `(departure, arrival)` cell. Each yielded `Result` is
/// independent — one input failing doesn't poison the rest.
///
/// For parallel evaluation, enable the `rayon` feature and use
/// [`lambert_par_iter`].
pub fn lambert_iter(
    inputs: &[LambertInput],
) -> impl Iterator<Item = Result<LambertSolutions, LambertError>> + '_ {
    inputs.iter().map(|input| input.solve())
}

/// Parallel batch iterator over Lambert inputs (Rayon-backed).
///
/// Same semantics as [`lambert_iter`], but evaluates inputs concurrently
/// across the Rayon thread pool. Available under the `rayon` feature.
#[cfg(feature = "rayon")]
#[allow(clippy::must_use_candidate)] // Rayon's ParallelIterator isn't #[must_use]; the caller is expected to chain `.for_each` / `.collect`.
pub fn lambert_par_iter(
    inputs: &[LambertInput],
) -> impl rayon::iter::ParallelIterator<Item = Result<LambertSolutions, LambertError>> + '_ {
    use rayon::prelude::*;
    inputs.par_iter().map(|input| input.solve())
}

/// Solve Lambert's boundary-value problem using Izzo's revisited algorithm.
///
/// Householder iteration over Lancaster's free parameter `x`, dispatching
/// across three TOF regimes (Battin / Lancaster–Blanchard / Lagrange) for
/// numerical stability. Mathematically frame-invariant — pass `r1` and
/// `r2` in any consistent inertial frame and the returned velocities are
/// in that same frame.
///
/// Returns the always-present single-revolution trajectory plus every
/// reachable multi-rev branch up to `revolutions.max()` (clamped at
/// [`MAX_MULTI_REV_PAIRS`]).
///
/// All quantities are dimensionally homogeneous — pass any consistent
/// unit system. The crate's docs and examples use km / km/s / s / km³/s².
///
/// # Invariants
///
/// All preconditions are validated at entry and returned as `Err(...)` on
/// violation — never panicked.
///
/// - `tof > 0`
/// - `mu > 0`
/// - `r1`, `r2`, `tof`, and `mu` are finite.
/// - `|r1| >= constants::MIN_POSITION_NORM`
/// - `|r2| >= constants::MIN_POSITION_NORM`
/// - Transfer angle ∉ {0, π}, equivalently
///   `|r1 × r2| / (|r1| · |r2|) >= constants::COLINEARITY_TOL`.
///
/// # Validity / near-degenerate behavior
///
/// - **Transfer angle near `0` or `π`** — the transfer plane is undefined;
///   returns [`LambertError::CollinearGeometry`]. Callers near these
///   boundaries should perturb one position by a tiny off-plane offset.
/// - **Near-parabolic (`|x − 1| ≤ 0.01`)** — the Lagrange and Lancaster TOF
///   formulations lose precision; the solver switches to Battin's
///   hypergeometric series (Izzo Eq. 20) automatically.
/// - **Hyperbolic transfers (`x > 1`)** — admitted on the single-rev branch;
///   multi-revolution solutions do not exist on a hyperbola and are silently
///   skipped.
/// - **Multi-rev infeasibility** — for `M ≥ 1`, the branch admits a solution
///   only when `tof ≥ T_min(M, λ)`. Higher-`M` branches are dropped from
///   the returned `multi` vector when their `T_min` exceeds the requested TOF.
///
/// # Errors
///
/// - [`LambertError::NonFiniteInput`] — any public scalar input or position
///   vector component is `NaN`, `+inf`, or `-inf`.
/// - [`LambertError::NonPositiveTimeOfFlight`] — `tof <= 0`.
/// - [`LambertError::NonPositiveMu`] — `mu <= 0`.
/// - [`LambertError::DegeneratePositionVector`] — `|r1|` or `|r2|` below
///   [`constants::MIN_POSITION_NORM`].
/// - [`LambertError::CollinearGeometry`] — `|r1 × r2| / (|r1| · |r2|)` below
///   [`constants::COLINEARITY_TOL`].
/// - [`LambertError::NoConvergence`] / [`LambertError::SingularDenominator`]
///   — Householder iteration failed.
pub fn lambert(
    r1: [f64; 3],
    r2: [f64; 3],
    tof: f64,
    mu: f64,
    way: TransferWay,
    revolutions: RevolutionBudget,
) -> Result<LambertSolutions, LambertError> {
    let geom = Geometry::from_inputs(r1, r2, tof, mu, way)?;
    let roots = find_xy(&geom, revolutions)?;
    Ok(reconstruct_solutions(&geom, &roots))
}

/// Like [`lambert`], but also returns the per-branch [`SolverDiagnostics`]
/// (iteration count and final Lancaster–Blanchard `x`).
///
/// The diagnostics structure mirrors the solutions structure 1:1 — `single`
/// matches `single`, `multi[i]` matches `multi[i]`.
///
/// # Invariants
///
/// Same as [`lambert`].
///
/// # Validity / near-degenerate behavior
///
/// Same as [`lambert`].
///
/// # Errors
///
/// Same as [`lambert`].
pub fn solve_with_diagnostics(
    r1: [f64; 3],
    r2: [f64; 3],
    tof: f64,
    mu: f64,
    way: TransferWay,
    revolutions: RevolutionBudget,
) -> Result<(LambertSolutions, LambertDiagnostics), LambertError> {
    let geom = Geometry::from_inputs(r1, r2, tof, mu, way)?;
    let roots = find_xy(&geom, revolutions)?;
    let solutions = reconstruct_solutions(&geom, &roots);
    let diagnostics = collect_diagnostics(&roots);
    Ok((solutions, diagnostics))
}

/// Solve both the short-way and long-way Lambert problems for the same
/// boundary inputs in one call.
///
/// Convenience for porkchop-plot and rendezvous-design callers that need
/// both traversal directions. Equivalent to calling [`lambert`] twice with
/// `TransferWay::Short` and `TransferWay::Long`; either half may have an
/// empty `multi` while the other is populated, since `T_min(M)` differs
/// between the two formulations.
///
/// # Invariants
///
/// Same as [`lambert`].
///
/// # Validity / near-degenerate behavior
///
/// Same as [`lambert`].
///
/// # Errors
///
/// Same as [`lambert`]. If either direction errors, the entire call errors —
/// the two ways share input validation.
pub fn lambert_both_ways(
    r1: [f64; 3],
    r2: [f64; 3],
    tof: f64,
    mu: f64,
    revolutions: RevolutionBudget,
) -> Result<BothWaysSolutions, LambertError> {
    let short = lambert(r1, r2, tof, mu, TransferWay::Short, revolutions)?;
    let long = lambert(r1, r2, tof, mu, TransferWay::Long, revolutions)?;
    Ok(BothWaysSolutions { short, long })
}

#[allow(clippy::similar_names)] // v_r1/v_r2/v_t1/v_t2 are radial/tangential velocity components at points 1/2 — Izzo §2.
fn reconstruct(geom: &Geometry, root: &Root) -> LambertSolution {
    // Velocity reconstruction (Izzo Algorithm 1).
    // Paper-named locals; `lambda*y ± x` are the elliptic-anomaly combos
    // that appear repeatedly in the radial/tangential decomposition.
    let lambda_y_minus_x = geom.lambda * root.y - root.x;
    let lambda_y_plus_x = geom.lambda * root.y + root.x;
    let tangential_num = geom.gamma * geom.sigma * (root.y + geom.lambda * root.x);

    let v_r1 = geom.gamma * (lambda_y_minus_x - geom.rho * lambda_y_plus_x) / geom.r1n;
    let v_r2 = -geom.gamma * (lambda_y_minus_x + geom.rho * lambda_y_plus_x) / geom.r2n;
    let v_t1 = tangential_num / geom.r1n;
    let v_t2 = tangential_num / geom.r2n;

    // v_radial · ir + v_tangential · it, component-wise.
    let v1 = [
        geom.ir1[0] * v_r1 + geom.it1[0] * v_t1,
        geom.ir1[1] * v_r1 + geom.it1[1] * v_t1,
        geom.ir1[2] * v_r1 + geom.it1[2] * v_t1,
    ];
    let v2 = [
        geom.ir2[0] * v_r2 + geom.it2[0] * v_t2,
        geom.ir2[1] * v_r2 + geom.it2[1] * v_t2,
        geom.ir2[2] * v_r2 + geom.it2[2] * v_t2,
    ];
    LambertSolution { v1, v2 }
}

fn reconstruct_solutions(geom: &Geometry, roots: &root_finding::Roots) -> LambertSolutions {
    let mut multi: ArrayVec<MultiRevPair, MAX_MULTI_REV_PAIRS> = ArrayVec::new();
    for pair in &roots.multi {
        // Capacity by construction — roots.multi is also bounded at MAX_MULTI_REV_PAIRS.
        let _ = multi.try_push(MultiRevPair {
            n_revs: pair.n_revs,
            long_period: reconstruct(geom, &pair.long_period),
            short_period: reconstruct(geom, &pair.short_period),
        });
    }
    LambertSolutions {
        single: reconstruct(geom, &roots.single),
        multi,
    }
}

fn diagnostics_of(root: &Root) -> SolverDiagnostics {
    SolverDiagnostics {
        iters: root.iters,
        lancaster_blanchard_x: root.x,
    }
}

fn collect_diagnostics(roots: &root_finding::Roots) -> LambertDiagnostics {
    let mut multi: ArrayVec<MultiRevPairDiagnostics, MAX_MULTI_REV_PAIRS> = ArrayVec::new();
    for pair in &roots.multi {
        let _ = multi.try_push(MultiRevPairDiagnostics {
            n_revs: pair.n_revs,
            long_period: diagnostics_of(&pair.long_period),
            short_period: diagnostics_of(&pair.short_period),
        });
    }
    LambertDiagnostics {
        single: diagnostics_of(&roots.single),
        multi,
    }
}

#[cfg(test)]
mod tests {
    #![allow(clippy::similar_names, clippy::unwrap_used)] // r/r1/r2 test scenario inputs follow paper convention; tests are exempt from the lib's no-unwrap rule.

    use super::*;
    use core::f64::consts::PI;
    use glam::DVec3;
    use lambert_izzo_test_support::bodies::{AU, MU_EARTH, MU_SUN};
    use lambert_izzo_test_support::kepler::propagate as kepler_propagate;
    use lambert_izzo_test_support::rand_unit_vec;

    fn approx(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    fn vec_sub_norm(a: [f64; 3], b: [f64; 3]) -> f64 {
        (DVec3::from_array(a) - DVec3::from_array(b)).length()
    }

    #[test]
    fn quarter_circle_leo() {
        // 90° transfer along a circular LEO at r = 7000 km.
        let r = 7000.0;
        let mu = MU_EARTH;
        let v_circ = (mu / r).sqrt();
        let period = 2.0 * PI * (r.powi(3) / mu).sqrt();

        let r1 = [r, 0.0, 0.0];
        let r2 = [0.0, r, 0.0];
        let sols = lambert(
            r1,
            r2,
            period / 4.0,
            mu,
            TransferWay::Short,
            RevolutionBudget::SingleOnly,
        )
        .unwrap();
        assert!(sols.multi.is_empty());
        assert!(vec_sub_norm(sols.single.v1, [0.0, v_circ, 0.0]) < 1e-9);
        assert!(vec_sub_norm(sols.single.v2, [-v_circ, 0.0, 0.0]) < 1e-9);
    }

    #[test]
    fn long_way_quarter_circle_leo() {
        // 270° transfer along the same circular LEO.
        let r = 7000.0;
        let mu = MU_EARTH;
        let v_circ = (mu / r).sqrt();
        let period = 2.0 * PI * (r.powi(3) / mu).sqrt();

        let r1 = [r, 0.0, 0.0];
        let r2 = [0.0, r, 0.0];
        let sols = lambert(
            r1,
            r2,
            3.0 * period / 4.0,
            mu,
            TransferWay::Long,
            RevolutionBudget::SingleOnly,
        )
        .unwrap();
        assert!(sols.multi.is_empty());
        assert!(vec_sub_norm(sols.single.v1, [0.0, -v_circ, 0.0]) < 1e-9);
        assert!(vec_sub_norm(sols.single.v2, [v_circ, 0.0, 0.0]) < 1e-9);
    }

    #[test]
    fn earth_mars_hohmann() {
        // Heliocentric Hohmann transfer Earth (1 AU) → Mars (1.524 AU).
        let mu = MU_SUN;
        let r1n = AU;
        let r2n = 1.524 * AU;
        let a = f64::midpoint(r1n, r2n);
        let tof = PI * (a.powi(3) / mu).sqrt();

        let r1 = [r1n, 0.0, 0.0];
        // 1 km off-plane to dodge the colinearity edge case.
        let r2 = [-r2n, 1.0, 0.0];
        let sols = lambert(
            r1,
            r2,
            tof,
            mu,
            TransferWay::Short,
            RevolutionBudget::SingleOnly,
        )
        .unwrap();
        // Periapsis velocity: vis-viva at r = 1 AU on the transfer ellipse.
        let v_peri = (mu * (2.0 / r1n - 1.0 / a)).sqrt();
        // Tolerance ~1 m/s — the off-plane perturbation moves things slightly.
        assert!(approx(DVec3::from_array(sols.single.v1).length(), v_peri, 1e-3));
    }

    #[test]
    fn multi_rev_branches() {
        // Long-tof Earth-orbit phasing — admits multiple multi-rev branches.
        let mu = MU_EARTH;
        let r1 = [8000.0, 0.0, 0.0];
        let r2 = [5600.0, 5600.0, 0.0];
        let period = 2.0 * PI * (8000.0_f64.powi(3) / mu).sqrt();
        let sols = lambert(
            r1,
            r2,
            5.0 * period,
            mu,
            TransferWay::Short,
            RevolutionBudget::up_to(3),
        )
        .unwrap();
        assert!(!sols.multi.is_empty(), "no multi-rev pairs returned");
        let r1n = DVec3::from_array(r1).length();
        for pair in &sols.multi {
            for s in [pair.long_period, pair.short_period] {
                let v1 = DVec3::from_array(s.v1);
                let energy = 0.5 * v1.dot(v1) - mu / r1n;
                assert!(energy.is_finite());
            }
        }
    }

    #[test]
    fn round_trip_kepler_check_single_rev() {
        // Propagate v1 with a universal-variable Kepler integrator and confirm
        // we land within numerical tolerance of r2.
        let mu = MU_EARTH;
        let r1 = [10_500.0, 1400.0, 700.0];
        let r2 = [-2800.0, 9100.0, -1400.0];
        let tof = 4500.0;
        let sols = lambert(
            r1,
            r2,
            tof,
            mu,
            TransferWay::Short,
            RevolutionBudget::SingleOnly,
        )
        .unwrap();
        let v1 = sols.single.v1;
        let r2_prop = kepler_propagate(r1, v1, tof, mu);
        let err = vec_sub_norm(r2_prop, r2);
        assert!(err < 1e-6, "kepler-roundtrip err = {err} km");
    }

    #[test]
    fn errors_on_non_positive_tof() {
        let r1 = [7000.0, 0.0, 0.0];
        let r2 = [0.0, 7000.0, 0.0];
        let err = lambert(
            r1,
            r2,
            0.0,
            MU_EARTH,
            TransferWay::Short,
            RevolutionBudget::SingleOnly,
        )
        .unwrap_err();
        assert!(matches!(err, LambertError::NonPositiveTimeOfFlight { tof } if tof == 0.0));
    }

    #[test]
    fn errors_on_zero_position_vector() {
        let r1 = [0.0, 0.0, 0.0];
        let r2 = [0.0, 7000.0, 0.0];
        let err = lambert(
            r1,
            r2,
            1000.0,
            MU_EARTH,
            TransferWay::Short,
            RevolutionBudget::SingleOnly,
        )
        .unwrap_err();
        assert!(matches!(
            err,
            LambertError::DegeneratePositionVector {
                position: Position::R1,
                ..
            }
        ));
    }

    #[test]
    fn errors_on_colinear_geometry() {
        let r1 = [7000.0, 0.0, 0.0];
        let r2 = [14_000.0, 0.0, 0.0];
        let err = lambert(
            r1,
            r2,
            1000.0,
            MU_EARTH,
            TransferWay::Short,
            RevolutionBudget::SingleOnly,
        )
        .unwrap_err();
        assert!(matches!(err, LambertError::CollinearGeometry { .. }));
    }

    #[test]
    fn errors_on_non_positive_mu() {
        let r1 = [7000.0, 0.0, 0.0];
        let r2 = [0.0, 7000.0, 0.0];
        let err = lambert(
            r1,
            r2,
            1000.0,
            0.0,
            TransferWay::Short,
            RevolutionBudget::SingleOnly,
        )
        .unwrap_err();
        assert!(matches!(err, LambertError::NonPositiveMu { mu } if mu == 0.0));
    }

    #[test]
    fn errors_on_non_finite_inputs() {
        let r1 = [7000.0, 0.0, 0.0];
        let r2 = [0.0, 7000.0, 0.0];

        let err = lambert(
            r1,
            r2,
            f64::NAN,
            MU_EARTH,
            TransferWay::Short,
            RevolutionBudget::SingleOnly,
        )
        .unwrap_err();
        assert!(matches!(
            err,
            LambertError::NonFiniteInput {
                parameter: NonFiniteParameter::Tof,
                ..
            }
        ));

        let err = lambert(
            [7000.0, f64::INFINITY, 0.0],
            r2,
            1000.0,
            MU_EARTH,
            TransferWay::Short,
            RevolutionBudget::SingleOnly,
        )
        .unwrap_err();
        assert!(matches!(
            err,
            LambertError::NonFiniteInput {
                parameter: NonFiniteParameter::R1Y,
                ..
            }
        ));
    }

    #[test]
    fn battin_regime_near_parabolic() {
        // GTO-like 90° transfer with TOF tuned so the converged x lands inside
        // the |x − 1| ≤ BATTIN_THRESHOLD (= 0.01) band, exercising the
        // hypergeometric series formulation in tof::x_to_tof_battin (Eq. 20).
        let mu = MU_EARTH;
        let r1 = [7000.0, 0.0, 0.0];
        let r2 = [0.0, 42_000.0, 0.0];
        let tof = 7200.0;
        let (sols, diag) = solve_with_diagnostics(
            r1,
            r2,
            tof,
            mu,
            TransferWay::Short,
            RevolutionBudget::SingleOnly,
        )
        .unwrap();
        assert!(
            (diag.single.lancaster_blanchard_x - 1.0).abs() < crate::constants::BATTIN_THRESHOLD,
            "x = {} not in Battin band [1 - {tol}, 1 + {tol}]",
            diag.single.lancaster_blanchard_x,
            tol = crate::constants::BATTIN_THRESHOLD,
        );
        let r2_prop = kepler_propagate(r1, sols.single.v1, tof, mu);
        let err = vec_sub_norm(r2_prop, r2);
        assert!(err < 1e-3, "Battin round-trip err = {err} km");
    }

    #[test]
    fn hyperbolic_transfer() {
        // Fast LEO → 200 000 km transfer; required v1 exceeds Earth escape,
        // landing in the hyperbolic branch (x > 1, positive specific energy).
        let mu = MU_EARTH;
        let r1 = [7000.0, 0.0, 0.0];
        let r2 = [0.0, 200_000.0, 0.0];
        let tof = 30_000.0;
        let (sols, diag) = solve_with_diagnostics(
            r1,
            r2,
            tof,
            mu,
            TransferWay::Short,
            RevolutionBudget::SingleOnly,
        )
        .unwrap();
        assert!(
            diag.single.lancaster_blanchard_x > 1.0,
            "expected hyperbolic (x > 1), got x = {}",
            diag.single.lancaster_blanchard_x
        );
        let v1 = DVec3::from_array(sols.single.v1);
        let energy = 0.5 * v1.dot(v1) - mu / DVec3::from_array(r1).length();
        assert!(
            energy > 0.0,
            "expected positive specific energy, got {energy}"
        );
        let r2_prop = kepler_propagate(r1, sols.single.v1, tof, mu);
        let err = vec_sub_norm(r2_prop, r2);
        assert!(err < 1e-3, "hyperbolic round-trip err = {err} km");
    }

    #[test]
    fn multi_rev_branches_distinct() {
        let mu = MU_EARTH;
        let r1 = [8000.0, 0.0, 0.0];
        let r2 = [5600.0, 5600.0, 0.0];
        let period = 2.0 * PI * (8000.0_f64.powi(3) / mu).sqrt();
        let tof = 5.0 * period;
        let (sols, diag) = solve_with_diagnostics(
            r1,
            r2,
            tof,
            mu,
            TransferWay::Short,
            RevolutionBudget::up_to(3),
        )
        .unwrap();

        assert!(!sols.multi.is_empty(), "no multi-rev pairs found");
        assert_eq!(sols.multi.len(), diag.multi.len());

        for (pair, diag_pair) in sols.multi.iter().zip(diag.multi.iter()) {
            assert_eq!(pair.n_revs, diag_pair.n_revs);
            assert!(
                diag_pair.long_period.lancaster_blanchard_x
                    < diag_pair.short_period.lancaster_blanchard_x,
                "M = {}: long-period x ({}) >= short-period x ({})",
                pair.n_revs,
                diag_pair.long_period.lancaster_blanchard_x,
                diag_pair.short_period.lancaster_blanchard_x,
            );
            for branch in [pair.long_period, pair.short_period] {
                let r2_prop = kepler_propagate(r1, branch.v1, tof, mu);
                let err = vec_sub_norm(r2_prop, r2);
                assert!(err < 1e-3, "M = {} branch round-trip err = {err} km", pair.n_revs);
            }
        }
    }

    #[test]
    fn solution_ordering_contract() {
        let mu = MU_EARTH;
        let r1 = [8000.0, 0.0, 0.0];
        let r2 = [5600.0, 5600.0, 0.0];
        let period = 2.0 * PI * (8000.0_f64.powi(3) / mu).sqrt();
        let tof = 5.0 * period;
        let sols = lambert(
            r1,
            r2,
            tof,
            mu,
            TransferWay::Short,
            RevolutionBudget::up_to(3),
        )
        .unwrap();

        let mut prev_m = 0_u32;
        for pair in &sols.multi {
            assert!(pair.n_revs > prev_m, "M strictly ascending across pairs");
            prev_m = pair.n_revs;
        }
    }

    #[cfg(feature = "serde")]
    #[test]
    fn serde_json_round_trip_preserves_solutions_and_errors() {
        // Solutions round-trip.
        let mu = MU_EARTH;
        let r1 = [8000.0, 0.0, 0.0];
        let r2 = [5600.0, 5600.0, 0.0];
        let period = 2.0 * PI * (8000.0_f64.powi(3) / mu).sqrt();
        let sols = lambert(
            r1,
            r2,
            5.0 * period,
            mu,
            TransferWay::Short,
            RevolutionBudget::up_to(2),
        )
        .unwrap();
        let json = serde_json::to_string(&sols).unwrap();
        let back: LambertSolutions = serde_json::from_str(&json).unwrap();
        assert_eq!(sols, back);
        assert!(!sols.multi.is_empty(), "test should exercise multi-rev branches");

        // Error round-trip — discriminated union via the default external tag.
        let err = LambertError::CollinearGeometry { sin_angle: 1e-20 };
        let err_json = serde_json::to_string(&err).unwrap();
        assert!(err_json.contains("CollinearGeometry"));
        let err_back: LambertError = serde_json::from_str(&err_json).unwrap();
        assert_eq!(err, err_back);
    }

    #[test]
    fn interop_with_nalgebra_and_glam_round_trips() {
        let mu = MU_EARTH;
        let r1_na = nalgebra::Vector3::new(7000.0, 0.0, 0.0);
        let r2_glam = glam::DVec3::new(0.0, 7000.0, 0.0);
        let r1: [f64; 3] = r1_na.into();
        let r2: [f64; 3] = r2_glam.to_array();
        let tof = PI / 2.0 * (7000.0_f64.powi(3) / mu).sqrt();
        let sols = lambert(
            r1,
            r2,
            tof,
            mu,
            TransferWay::Short,
            RevolutionBudget::SingleOnly,
        )
        .unwrap();
        let v1_back_na: nalgebra::Vector3<f64> = sols.single.v1.into();
        let v2_back_glam = glam::DVec3::from_array(sols.single.v2);
        assert!(v1_back_na.iter().all(|c| c.is_finite()));
        assert!(v2_back_glam.to_array().iter().all(|c| c.is_finite()));
    }

    #[test]
    fn both_ways_returns_independent_halves() {
        let mu = MU_EARTH;
        let r: f64 = 7000.0;
        let r1 = [r, 0.0, 0.0];
        let r2 = [0.0, r, 0.0];
        let period = 2.0 * PI * (r.powi(3) / mu).sqrt();
        let tof = period / 4.0;

        let both = lambert_both_ways(r1, r2, tof, mu, RevolutionBudget::SingleOnly).unwrap();
        let short = lambert(
            r1,
            r2,
            tof,
            mu,
            TransferWay::Short,
            RevolutionBudget::SingleOnly,
        )
        .unwrap();
        let long = lambert(
            r1,
            r2,
            3.0 * period / 4.0,
            mu,
            TransferWay::Long,
            RevolutionBudget::SingleOnly,
        )
        .unwrap();

        assert_eq!(both.short, short);
        assert!(both.long.single.v1.iter().all(|c| c.is_finite()));
        let v_circ = (mu / r).sqrt();
        assert!(vec_sub_norm(long.single.v1, [0.0, -v_circ, 0.0]) < 1e-9);
    }

    #[test]
    fn kepler_roundtrip_random_single_rev() {
        use rand::{Rng, SeedableRng};
        use rand_chacha::ChaCha20Rng;
        use rand_distr::Uniform;

        let mu = MU_EARTH;
        let mut rng = ChaCha20Rng::seed_from_u64(0xC0FF_EE42);
        let radius = Uniform::new(3500.0, 28_000.0);
        let tof_dist = Uniform::new(100.0, 50_000.0);

        let mut max_rel_err = 0.0_f64;
        let mut good_count = 0_u32;
        let mut lambert_ok = 0_u32;
        for _ in 0..1000 {
            let (r1_mag, r2_mag) = (rng.sample(radius), rng.sample(radius));
            let r1 = (DVec3::from_array(rand_unit_vec(&mut rng)) * r1_mag).to_array();
            let r2 = (DVec3::from_array(rand_unit_vec(&mut rng)) * r2_mag).to_array();
            let tof = rng.sample(tof_dist);
            let way = if rng.gen_bool(0.5) {
                TransferWay::Long
            } else {
                TransferWay::Short
            };
            let Ok(sols) = lambert(r1, r2, tof, mu, way, RevolutionBudget::SingleOnly) else {
                continue;
            };
            lambert_ok += 1;
            let r2_prop = kepler_propagate(r1, sols.single.v1, tof, mu);
            let rel = vec_sub_norm(r2_prop, r2) / DVec3::from_array(r2).length();
            if rel.is_finite() {
                max_rel_err = max_rel_err.max(rel);
                good_count += 1;
            }
        }
        assert!(lambert_ok > 950, "too many Lambert failures: {lambert_ok}/1000");
        assert!(good_count > 500, "too few converged round-trips: {good_count}/{lambert_ok}");
        assert!(max_rel_err < 1e-6, "max relative round-trip err = {max_rel_err:.3e} over {good_count} trials");
    }

    #[test]
    fn kepler_roundtrip_random_multi_rev() {
        use rand::{Rng, SeedableRng};
        use rand_chacha::ChaCha20Rng;
        use rand_distr::Uniform;

        let mu = MU_EARTH;
        let mut rng = ChaCha20Rng::seed_from_u64(0xBEEF_DEAD);
        let radius = Uniform::new(5600.0, 10_500.0);
        let tof_dist = Uniform::new(10_000.0, 250_000.0);

        let mut max_rel_err = 0.0_f64;
        let mut branches = 0_u32;
        let mut good_count = 0_u32;
        for _ in 0..500 {
            let (r1_mag, r2_mag) = (rng.sample(radius), rng.sample(radius));
            let r1 = (DVec3::from_array(rand_unit_vec(&mut rng)) * r1_mag).to_array();
            let r2 = (DVec3::from_array(rand_unit_vec(&mut rng)) * r2_mag).to_array();
            let tof = rng.sample(tof_dist);
            let Ok(sols) = lambert(r1, r2, tof, mu, TransferWay::Short, RevolutionBudget::up_to(3)) else {
                continue;
            };
            let mut iter_branch = |s: LambertSolution| {
                let r2_prop = kepler_propagate(r1, s.v1, tof, mu);
                let rel = vec_sub_norm(r2_prop, r2) / DVec3::from_array(r2).length();
                if rel.is_finite() {
                    max_rel_err = max_rel_err.max(rel);
                    good_count += 1;
                }
                branches += 1;
            };
            iter_branch(sols.single);
            for pair in &sols.multi {
                iter_branch(pair.long_period);
                iter_branch(pair.short_period);
            }
        }
        assert!(branches > 500, "expected branches > 500, got {branches}");
        assert!(good_count > 300, "too few converged round-trips: {good_count}/{branches}");
        assert!(max_rel_err < 1e-5, "max relative round-trip err = {max_rel_err:.3e} over {good_count} branches");
    }
}
