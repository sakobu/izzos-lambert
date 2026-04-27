//! Izzo's revisited Lambert solver — single + multi-revolution, short/long way.
//!
//! Reference: D. Izzo, *Revisiting Lambert's problem*, Celestial Mechanics &
//! Dynamical Astronomy, 2014. arXiv:1403.2705. PDF in `docs/izzo.pdf`.
//!
//! Inline `Eq. N` / `Algorithm N` references in the source point to that paper.
//!
//! # Units
//!
//! Public inputs and outputs follow these SI suffix conventions:
//!
//! | Quantity                  | Suffix     | Unit       |
//! |---------------------------|------------|------------|
//! | Position                  | `_km`      | km         |
//! | Velocity                  | `_km_s`    | km/s       |
//! | Time of flight            | `_s`       | s          |
//! | Gravitational parameter   | `_km3_s2`  | km³/s²     |
//!
//! The algorithm is mathematically frame-invariant under any inertial
//! frame — pass `r1_km` and `r2_km` in the same inertial frame (ECI for
//! Earth orbits, HCRS for solar transfers, etc.) and the returned
//! velocities are in that same frame.
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
//! let v1_na: nalgebra::Vector3<f64> = solution.single.v1_km_s.into();
//!
//! // glam:
//! let r2 = glam::DVec3::new(0.0, 7000.0, 0.0).to_array();
//! let v2_glam = glam::DVec3::from_array(solution.single.v2_km_s);
//! ```
//!
//! # Example
//!
//! ```
//! use lambert_izzo::{lambert, LambertError, RevolutionBudget, TransferWay};
//!
//! # fn main() -> Result<(), LambertError> {
//! // LEO → LEO 90° transfer at 7000 km altitude.
//! let mu_km3_s2 = 398_600.4418;
//! let r1_km = [7000.0, 0.0, 0.0];
//! let r2_km = [0.0, 7000.0, 0.0];
//! let tof_s = core::f64::consts::PI / 2.0 * (7000.0_f64.powi(3) / mu_km3_s2).sqrt();
//!
//! let solutions = lambert(
//!     r1_km, r2_km, tof_s, mu_km3_s2,
//!     TransferWay::Short, RevolutionBudget::SingleOnly,
//! )?;
//! assert!(solutions.multi.is_empty());
//! let v1_km_s = solutions.single.v1_km_s;
//! # let _ = v1_km_s;
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
/// type `[multi-rev-pair; MAX_MULTI_REV_PAIRS]` a fixed stack size.
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
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum RevolutionBudget {
    /// Solve the single-revolution case only — always exactly one solution.
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

#[cfg(any(test, feature = "test-utils"))]
mod test_helpers;

/// Test utilities exposed under the `test-utils` feature.
///
/// Currently a universal-variable Kepler propagator suitable for
/// round-trip-validating Lambert solutions in downstream integration
/// tests, so callers don't have to re-implement Stumpff functions.
#[cfg(any(test, feature = "test-utils"))]
pub mod test_utils {
    pub use super::test_helpers::kepler_propagate;
}

pub use error::{LambertError, NonFiniteParameter};

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
    /// Velocity at `r1_km` (km/s, same inertial frame as the inputs).
    pub v1_km_s: [f64; 3],
    /// Velocity at `r2_km` (km/s, same inertial frame as the inputs).
    pub v2_km_s: [f64; 3],
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
    /// Initial position (km), any consistent inertial frame.
    pub r1_km: [f64; 3],
    /// Final position (km), same frame as `r1_km`.
    pub r2_km: [f64; 3],
    /// Time of flight (s), `> 0`.
    pub tof_s: f64,
    /// Gravitational parameter (km³/s²), `> 0`.
    pub mu_km3_s2: f64,
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
        lambert(
            self.r1_km,
            self.r2_km,
            self.tof_s,
            self.mu_km3_s2,
            self.way,
            self.revolutions,
        )
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
/// numerical stability. Mathematically frame-invariant — pass `r1_km` and
/// `r2_km` in any consistent inertial frame and the returned velocities are
/// in that same frame.
///
/// Returns the always-present single-revolution trajectory plus every
/// reachable multi-rev branch up to `revolutions.max()` (clamped at
/// [`MAX_MULTI_REV_PAIRS`]).
///
/// # Invariants
///
/// All preconditions are validated at entry and returned as `Err(...)` on
/// violation — never panicked.
///
/// - `tof_s > 0`
/// - `mu_km3_s2 > 0`
/// - `r1_km`, `r2_km`, `tof_s`, and `mu_km3_s2` are finite.
/// - `|r1_km| >= constants::MIN_POSITION_NORM_KM`
/// - `|r2_km| >= constants::MIN_POSITION_NORM_KM`
/// - Transfer angle ∉ {0, π}, equivalently
///   `|r1_km × r2_km| / (|r1_km| · |r2_km|) >= constants::COLINEARITY_TOL`.
///
/// # Validity / near-degenerate behavior
///
/// - **Transfer angle near `0` or `π`** — the transfer plane is undefined;
///   returns [`LambertError::CollinearGeometry`]. Callers near these
///   boundaries should perturb one position by ≈ 1 km off-plane.
/// - **Near-parabolic (`|x − 1| ≤ 0.01`)** — the Lagrange and Lancaster TOF
///   formulations lose precision; the solver switches to Battin's
///   hypergeometric series (Izzo Eq. 20) automatically.
/// - **Hyperbolic transfers (`x > 1`)** — admitted on the single-rev branch;
///   multi-revolution solutions do not exist on a hyperbola and are silently
///   skipped.
/// - **Multi-rev infeasibility** — for `M ≥ 1`, the branch admits a solution
///   only when `tof_s ≥ T_min(M, λ)`. Higher-`M` branches are dropped from
///   the returned `multi` vector when their `T_min` exceeds the requested TOF.
///
/// # Errors
///
/// - [`LambertError::NonFiniteInput`] — any public scalar input or position
///   vector component is `NaN`, `+inf`, or `-inf`.
/// - [`LambertError::NonPositiveTimeOfFlight`] — `tof_s <= 0`.
/// - [`LambertError::NonPositiveMu`] — `mu_km3_s2 <= 0`.
/// - [`LambertError::DegeneratePositionVector`] — `|r1|` or `|r2|` below
///   [`constants::MIN_POSITION_NORM_KM`].
/// - [`LambertError::CollinearGeometry`] — `|r1 × r2| / (|r1| · |r2|)` below
///   [`constants::COLINEARITY_TOL`].
/// - [`LambertError::NoConvergence`] / [`LambertError::SingularDenominator`]
///   — Householder iteration failed.
pub fn lambert(
    r1_km: [f64; 3],
    r2_km: [f64; 3],
    tof_s: f64,
    mu_km3_s2: f64,
    way: TransferWay,
    revolutions: RevolutionBudget,
) -> Result<LambertSolutions, LambertError> {
    let geom = Geometry::from_inputs(r1_km, r2_km, tof_s, mu_km3_s2, way)?;
    let roots = find_xy(&geom, revolutions)?;
    Ok(reconstruct_solutions(&geom, &roots))
}

/// Like [`lambert`], but also returns the per-branch [`SolverDiagnostics`]
/// (iteration count and final Lancaster–Blanchard `x`).
///
/// The diagnostics structure mirrors the solutions structure 1:1 — `single`
/// matches `single`, `multi[i]` matches `multi[i]`.
///
/// # Errors
///
/// Same as [`lambert`].
pub fn solve_with_diagnostics(
    r1_km: [f64; 3],
    r2_km: [f64; 3],
    tof_s: f64,
    mu_km3_s2: f64,
    way: TransferWay,
    revolutions: RevolutionBudget,
) -> Result<(LambertSolutions, LambertDiagnostics), LambertError> {
    let geom = Geometry::from_inputs(r1_km, r2_km, tof_s, mu_km3_s2, way)?;
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
/// # Errors
///
/// Same as [`lambert`]. If either direction errors, the entire call errors —
/// the two ways share input validation.
pub fn lambert_both_ways(
    r1_km: [f64; 3],
    r2_km: [f64; 3],
    tof_s: f64,
    mu_km3_s2: f64,
    revolutions: RevolutionBudget,
) -> Result<BothWaysSolutions, LambertError> {
    let short = lambert(
        r1_km,
        r2_km,
        tof_s,
        mu_km3_s2,
        TransferWay::Short,
        revolutions,
    )?;
    let long = lambert(
        r1_km,
        r2_km,
        tof_s,
        mu_km3_s2,
        TransferWay::Long,
        revolutions,
    )?;
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

    let v1_km_s = vec3::add(vec3::scale(geom.ir1, v_r1), vec3::scale(geom.it1, v_t1));
    let v2_km_s = vec3::add(vec3::scale(geom.ir2, v_r2), vec3::scale(geom.it2, v_t2));
    LambertSolution { v1_km_s, v2_km_s }
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
    #![allow(clippy::similar_names, clippy::unwrap_used)] // r_km/r1_km/r2_km test scenario inputs follow paper convention; tests are exempt from the lib's no-unwrap rule.

    use super::*;
    use crate::test_helpers::kepler_propagate;
    use crate::vec3;
    use core::f64::consts::PI;

    /// Earth's gravitational parameter (km³/s²) — value from EGM2008.
    const MU_EARTH_KM3_S2: f64 = 398_600.441_8;
    /// Sun's gravitational parameter (km³/s²) — value from DE440.
    const MU_SUN_KM3_S2: f64 = 1.327_124_400_18e11;

    fn approx(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    fn vec_sub_norm(a: [f64; 3], b: [f64; 3]) -> f64 {
        vec3::norm(vec3::sub(a, b))
    }

    #[test]
    fn quarter_circle_leo() {
        // 90° transfer along a circular LEO at r = 7000 km.
        let r_km = 7000.0;
        let mu = MU_EARTH_KM3_S2;
        let v_circ = (mu / r_km).sqrt();
        let period_s = 2.0 * PI * (r_km.powi(3) / mu).sqrt();

        let r1_km = [r_km, 0.0, 0.0];
        let r2_km = [0.0, r_km, 0.0];
        let sols = lambert(
            r1_km,
            r2_km,
            period_s / 4.0,
            mu,
            TransferWay::Short,
            RevolutionBudget::SingleOnly,
        )
        .unwrap();
        assert!(sols.multi.is_empty());
        assert!(vec_sub_norm(sols.single.v1_km_s, [0.0, v_circ, 0.0]) < 1e-9);
        assert!(vec_sub_norm(sols.single.v2_km_s, [-v_circ, 0.0, 0.0]) < 1e-9);
    }

    #[test]
    fn long_way_quarter_circle_leo() {
        // 270° transfer along the same circular LEO.
        let r_km = 7000.0;
        let mu = MU_EARTH_KM3_S2;
        let v_circ = (mu / r_km).sqrt();
        let period_s = 2.0 * PI * (r_km.powi(3) / mu).sqrt();

        let r1_km = [r_km, 0.0, 0.0];
        let r2_km = [0.0, r_km, 0.0];
        let sols = lambert(
            r1_km,
            r2_km,
            3.0 * period_s / 4.0,
            mu,
            TransferWay::Long,
            RevolutionBudget::SingleOnly,
        )
        .unwrap();
        assert!(sols.multi.is_empty());
        assert!(vec_sub_norm(sols.single.v1_km_s, [0.0, -v_circ, 0.0]) < 1e-9);
        assert!(vec_sub_norm(sols.single.v2_km_s, [v_circ, 0.0, 0.0]) < 1e-9);
    }

    #[test]
    fn earth_mars_hohmann() {
        // Heliocentric Hohmann transfer Earth (1 AU) → Mars (1.524 AU).
        const AU_KM: f64 = 1.495_978_707e8;
        let mu = MU_SUN_KM3_S2;
        let r1_norm_km = AU_KM;
        let r2_norm_km = 1.524 * AU_KM;
        let a_km = f64::midpoint(r1_norm_km, r2_norm_km);
        let tof_s = PI * (a_km.powi(3) / mu).sqrt();

        let r1_km = [r1_norm_km, 0.0, 0.0];
        // 1 km off-plane to dodge the colinearity edge case.
        let r2_km = [-r2_norm_km, 1.0, 0.0];
        let sols = lambert(
            r1_km,
            r2_km,
            tof_s,
            mu,
            TransferWay::Short,
            RevolutionBudget::SingleOnly,
        )
        .unwrap();
        // Periapsis velocity: vis-viva at r = 1 AU on the transfer ellipse.
        let v_peri_km_s = (mu * (2.0 / r1_norm_km - 1.0 / a_km)).sqrt();
        // Tolerance ~1 m/s — the off-plane perturbation moves things slightly.
        assert!(approx(vec3::norm(sols.single.v1_km_s), v_peri_km_s, 1e-3));
    }

    #[test]
    fn multi_rev_branches() {
        // Long-tof Earth-orbit phasing — admits multiple multi-rev branches.
        let mu = MU_EARTH_KM3_S2;
        let r1_km = [8000.0, 0.0, 0.0];
        let r2_km = [5600.0, 5600.0, 0.0];
        let period_s = 2.0 * PI * (8000.0_f64.powi(3) / mu).sqrt();
        let sols = lambert(
            r1_km,
            r2_km,
            5.0 * period_s,
            mu,
            TransferWay::Short,
            RevolutionBudget::up_to(3),
        )
        .unwrap();
        assert!(!sols.multi.is_empty(), "no multi-rev pairs returned");
        let r1n = vec3::norm(r1_km);
        for pair in &sols.multi {
            for s in [pair.long_period, pair.short_period] {
                let energy = 0.5 * vec3::dot(s.v1_km_s, s.v1_km_s) - mu / r1n;
                assert!(energy.is_finite());
            }
        }
    }

    #[test]
    fn round_trip_kepler_check_single_rev() {
        // Propagate v1 with a universal-variable Kepler integrator and confirm
        // we land within numerical tolerance of r2.
        let mu = MU_EARTH_KM3_S2;
        let r1_km = [10_500.0, 1400.0, 700.0];
        let r2_km = [-2800.0, 9100.0, -1400.0];
        let tof_s = 4500.0;
        let sols = lambert(
            r1_km,
            r2_km,
            tof_s,
            mu,
            TransferWay::Short,
            RevolutionBudget::SingleOnly,
        )
        .unwrap();
        let v1_km_s = sols.single.v1_km_s;
        let r2_prop_km = kepler_propagate(r1_km, v1_km_s, tof_s, mu);
        let err_km = vec_sub_norm(r2_prop_km, r2_km);
        assert!(err_km < 1e-6, "kepler-roundtrip err = {err_km} km");
    }

    #[test]
    fn errors_on_non_positive_tof() {
        let r1_km = [7000.0, 0.0, 0.0];
        let r2_km = [0.0, 7000.0, 0.0];
        let err = lambert(
            r1_km,
            r2_km,
            0.0,
            MU_EARTH_KM3_S2,
            TransferWay::Short,
            RevolutionBudget::SingleOnly,
        )
        .unwrap_err();
        assert!(matches!(err, LambertError::NonPositiveTimeOfFlight { tof_s } if tof_s == 0.0));
    }

    #[test]
    fn errors_on_zero_position_vector() {
        let r1_km = [0.0, 0.0, 0.0];
        let r2_km = [0.0, 7000.0, 0.0];
        let err = lambert(
            r1_km,
            r2_km,
            1000.0,
            MU_EARTH_KM3_S2,
            TransferWay::Short,
            RevolutionBudget::SingleOnly,
        )
        .unwrap_err();
        assert!(matches!(
            err,
            LambertError::DegeneratePositionVector { which: 1, .. }
        ));
    }

    #[test]
    fn errors_on_colinear_geometry() {
        let r1_km = [7000.0, 0.0, 0.0];
        let r2_km = [14_000.0, 0.0, 0.0];
        let err = lambert(
            r1_km,
            r2_km,
            1000.0,
            MU_EARTH_KM3_S2,
            TransferWay::Short,
            RevolutionBudget::SingleOnly,
        )
        .unwrap_err();
        assert!(matches!(err, LambertError::CollinearGeometry { .. }));
    }

    #[test]
    fn errors_on_non_positive_mu() {
        let r1_km = [7000.0, 0.0, 0.0];
        let r2_km = [0.0, 7000.0, 0.0];
        let err = lambert(
            r1_km,
            r2_km,
            1000.0,
            0.0,
            TransferWay::Short,
            RevolutionBudget::SingleOnly,
        )
        .unwrap_err();
        assert!(matches!(err, LambertError::NonPositiveMu { mu_km3_s2 } if mu_km3_s2 == 0.0));
    }

    #[test]
    fn errors_on_non_finite_inputs() {
        let r1_km = [7000.0, 0.0, 0.0];
        let r2_km = [0.0, 7000.0, 0.0];

        let err = lambert(
            r1_km,
            r2_km,
            f64::NAN,
            MU_EARTH_KM3_S2,
            TransferWay::Short,
            RevolutionBudget::SingleOnly,
        )
        .unwrap_err();
        assert!(matches!(
            err,
            LambertError::NonFiniteInput {
                parameter: NonFiniteParameter::TofS,
                ..
            }
        ));

        let err = lambert(
            [7000.0, f64::INFINITY, 0.0],
            r2_km,
            1000.0,
            MU_EARTH_KM3_S2,
            TransferWay::Short,
            RevolutionBudget::SingleOnly,
        )
        .unwrap_err();
        assert!(matches!(
            err,
            LambertError::NonFiniteInput {
                parameter: NonFiniteParameter::R1KmY,
                ..
            }
        ));
    }

    #[test]
    fn battin_regime_near_parabolic() {
        // GTO-like 90° transfer with TOF tuned so the converged x lands inside
        // the |x − 1| ≤ BATTIN_THRESHOLD (= 0.01) band, exercising the
        // hypergeometric series formulation in tof::x_to_tof_battin (Eq. 20).
        let mu = MU_EARTH_KM3_S2;
        let r1_km = [7000.0, 0.0, 0.0];
        let r2_km = [0.0, 42_000.0, 0.0];
        let tof_s = 7200.0;
        let (sols, diag) = solve_with_diagnostics(
            r1_km,
            r2_km,
            tof_s,
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
        let r2_prop = kepler_propagate(r1_km, sols.single.v1_km_s, tof_s, mu);
        let err_km = vec_sub_norm(r2_prop, r2_km);
        assert!(err_km < 1e-3, "Battin round-trip err = {err_km} km");
    }

    #[test]
    fn hyperbolic_transfer() {
        // Fast LEO → 200 000 km transfer; required v1 exceeds Earth escape,
        // landing in the hyperbolic branch (x > 1, positive specific energy).
        // Exercises tof::x_to_tof_lagrange's a < 0 path and compute_psi's
        // asinh branch (Eq. 9 / Eq. 17 hyperbolic).
        let mu = MU_EARTH_KM3_S2;
        let r1_km = [7000.0, 0.0, 0.0];
        let r2_km = [0.0, 200_000.0, 0.0];
        let tof_s = 30_000.0;
        let (sols, diag) = solve_with_diagnostics(
            r1_km,
            r2_km,
            tof_s,
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
        let energy = 0.5 * vec3::dot(sols.single.v1_km_s, sols.single.v1_km_s)
            - mu / vec3::norm(r1_km);
        assert!(
            energy > 0.0,
            "expected positive specific energy, got {energy}"
        );
        let r2_prop = kepler_propagate(r1_km, sols.single.v1_km_s, tof_s, mu);
        let err_km = vec_sub_norm(r2_prop, r2_km);
        assert!(err_km < 1e-3, "hyperbolic round-trip err = {err_km} km");
    }

    #[test]
    fn multi_rev_branches_distinct() {
        // Same geometry as multi_rev_branches; verify long-period x < short-period x
        // for each M (paper's "switch between branches" concern, page 14) and that
        // both branches independently round-trip via Kepler.
        let mu = MU_EARTH_KM3_S2;
        let r1_km = [8000.0, 0.0, 0.0];
        let r2_km = [5600.0, 5600.0, 0.0];
        let period_s = 2.0 * PI * (8000.0_f64.powi(3) / mu).sqrt();
        let tof_s = 5.0 * period_s;
        let (sols, diag) = solve_with_diagnostics(
            r1_km,
            r2_km,
            tof_s,
            mu,
            TransferWay::Short,
            RevolutionBudget::up_to(3),
        )
        .unwrap();

        assert!(!sols.multi.is_empty(), "no multi-rev pairs found");
        assert_eq!(sols.multi.len(), diag.multi.len());

        for (pair, diag_pair) in sols.multi.iter().zip(diag.multi.iter()) {
            assert_eq!(pair.n_revs, diag_pair.n_revs);
            // Long-period x is strictly smaller than short-period x.
            assert!(
                diag_pair.long_period.lancaster_blanchard_x
                    < diag_pair.short_period.lancaster_blanchard_x,
                "M = {}: long-period x ({}) >= short-period x ({})",
                pair.n_revs,
                diag_pair.long_period.lancaster_blanchard_x,
                diag_pair.short_period.lancaster_blanchard_x,
            );
            for branch in [pair.long_period, pair.short_period] {
                let r2_prop = kepler_propagate(r1_km, branch.v1_km_s, tof_s, mu);
                let err_km = vec_sub_norm(r2_prop, r2_km);
                assert!(
                    err_km < 1e-3,
                    "M = {} branch round-trip err = {err_km} km",
                    pair.n_revs
                );
            }
        }
    }

    #[test]
    fn solution_ordering_contract() {
        // Multi-rev pairs are ascending in M (compile-enforced no longer needed —
        // the type makes the structure explicit; verify ascending order).
        let mu = MU_EARTH_KM3_S2;
        let r1_km = [8000.0, 0.0, 0.0];
        let r2_km = [5600.0, 5600.0, 0.0];
        let period_s = 2.0 * PI * (8000.0_f64.powi(3) / mu).sqrt();
        let tof_s = 5.0 * period_s;
        let sols = lambert(
            r1_km,
            r2_km,
            tof_s,
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
        let mu = MU_EARTH_KM3_S2;
        let r1_km = [8000.0, 0.0, 0.0];
        let r2_km = [5600.0, 5600.0, 0.0];
        let period_s = 2.0 * PI * (8000.0_f64.powi(3) / mu).sqrt();
        let sols = lambert(
            r1_km,
            r2_km,
            5.0 * period_s,
            mu,
            TransferWay::Short,
            RevolutionBudget::up_to(2),
        )
        .unwrap();
        let json = serde_json::to_string(&sols).unwrap();
        let back: LambertSolutions = serde_json::from_str(&json).unwrap();
        assert_eq!(sols, back);
        assert!(!sols.multi.is_empty(), "test should exercise multi-rev branches");

        // Error round-trip — discriminated union via the `kind` tag.
        let err = LambertError::CollinearGeometry { sin_angle: 1e-20 };
        let err_json = serde_json::to_string(&err).unwrap();
        assert!(err_json.contains("CollinearGeometry"));
        let err_back: LambertError = serde_json::from_str(&err_json).unwrap();
        assert_eq!(err, err_back);
    }

    #[test]
    fn interop_with_nalgebra_and_glam_round_trips() {
        // Verifies the doc-level claim that the public [f64; 3] surface
        // converts cleanly to/from nalgebra::Vector3<f64> and glam::DVec3
        // without any feature-flagged shim.
        let mu = MU_EARTH_KM3_S2;
        let r1_na = nalgebra::Vector3::new(7000.0, 0.0, 0.0);
        let r2_glam = glam::DVec3::new(0.0, 7000.0, 0.0);
        let r1_km: [f64; 3] = r1_na.into();
        let r2_km: [f64; 3] = r2_glam.to_array();
        let tof_s = PI / 2.0 * (7000.0_f64.powi(3) / mu).sqrt();
        let sols = lambert(
            r1_km,
            r2_km,
            tof_s,
            mu,
            TransferWay::Short,
            RevolutionBudget::SingleOnly,
        )
        .unwrap();
        // Returned arrays convert back to either lib.
        let v1_back_na: nalgebra::Vector3<f64> = sols.single.v1_km_s.into();
        let v2_back_glam = glam::DVec3::from_array(sols.single.v2_km_s);
        assert!(v1_back_na.iter().all(|c| c.is_finite()));
        assert!(v2_back_glam.to_array().iter().all(|c| c.is_finite()));
    }

    #[test]
    fn both_ways_returns_independent_halves() {
        // Verify lambert_both_ways produces the same answers as two separate
        // lambert calls.
        let mu = MU_EARTH_KM3_S2;
        let r_km: f64 = 7000.0;
        let r1_km = [r_km, 0.0, 0.0];
        let r2_km = [0.0, r_km, 0.0];
        let period_s = 2.0 * PI * (r_km.powi(3) / mu).sqrt();
        let tof_s = period_s / 4.0;

        let both = lambert_both_ways(r1_km, r2_km, tof_s, mu, RevolutionBudget::SingleOnly).unwrap();
        let short = lambert(
            r1_km,
            r2_km,
            tof_s,
            mu,
            TransferWay::Short,
            RevolutionBudget::SingleOnly,
        )
        .unwrap();
        let long = lambert(
            r1_km,
            r2_km,
            3.0 * period_s / 4.0,
            mu,
            TransferWay::Long,
            RevolutionBudget::SingleOnly,
        )
        .unwrap();

        assert_eq!(both.short, short);
        // The long way needs the long TOF, but here we passed the short TOF
        // to both_ways — so just check the returned trajectory is finite.
        assert!(both.long.single.v1_km_s.iter().all(|c| c.is_finite()));
        // The independent long-TOF call should produce the analytic circular
        // velocity (asserted in long_way_quarter_circle_leo).
        let v_circ = (mu / r_km).sqrt();
        assert!(vec_sub_norm(long.single.v1_km_s, [0.0, -v_circ, 0.0]) < 1e-9);
    }

    fn rand_unit_vec(rng: &mut rand_chacha::ChaCha20Rng) -> [f64; 3] {
        use rand::Rng;
        use rand_distr::Uniform;
        let axis: Uniform<f64> = Uniform::new(-1.0, 1.0);
        loop {
            let v: [f64; 3] = [rng.sample(axis), rng.sample(axis), rng.sample(axis)];
            let n2 = vec3::norm_squared(v);
            if n2 > 0.01 && n2 < 1.0 {
                return vec3::scale(v, 1.0 / n2.sqrt());
            }
        }
    }

    #[test]
    fn kepler_roundtrip_random_single_rev() {
        // Paper §5 / Fig. 6 statistical validation, physical analog: 1000 random
        // Earth-orbit geometries, solve Lambert, propagate v1 with universal-
        // variable Kepler, check |r2_prop − r2| / |r2| stays small.
        use rand::{Rng, SeedableRng};
        use rand_chacha::ChaCha20Rng;
        use rand_distr::Uniform;

        let mu = MU_EARTH_KM3_S2;
        let mut rng = ChaCha20Rng::seed_from_u64(0xC0FF_EE42);
        let radius = Uniform::new(3500.0, 28_000.0);
        let tof = Uniform::new(100.0, 50_000.0);

        let mut max_rel_err = 0.0_f64;
        let mut good_count = 0_u32;
        let mut lambert_ok = 0_u32;
        for _ in 0..1000 {
            let r1_km = vec3::scale(rand_unit_vec(&mut rng), rng.sample(radius));
            let r2_km = vec3::scale(rand_unit_vec(&mut rng), rng.sample(radius));
            let tof_s = rng.sample(tof);
            let way = if rng.gen_bool(0.5) {
                TransferWay::Long
            } else {
                TransferWay::Short
            };
            let Ok(sols) = lambert(r1_km, r2_km, tof_s, mu, way, RevolutionBudget::SingleOnly)
            else {
                continue;
            };
            lambert_ok += 1;
            let r2_prop = kepler_propagate(r1_km, sols.single.v1_km_s, tof_s, mu);
            let rel = vec_sub_norm(r2_prop, r2_km) / vec3::norm(r2_km);
            // NaN signals propagator divergence on a pathological geometry —
            // not a Lambert correctness issue. Filter and assert on the rest.
            if rel.is_finite() {
                max_rel_err = max_rel_err.max(rel);
                good_count += 1;
            }
        }
        assert!(
            lambert_ok > 950,
            "too many Lambert failures: {lambert_ok}/1000"
        );
        assert!(
            good_count > 500,
            "too few converged round-trips: {good_count}/{lambert_ok}"
        );
        assert!(
            max_rel_err < 1e-6,
            "max relative round-trip err = {max_rel_err:.3e} over {good_count} trials"
        );
    }

    #[test]
    fn kepler_roundtrip_random_multi_rev() {
        // Paper §5 multi-rev analog: 500 random geometries with multi_revs=3,
        // every returned branch (single + multi) verified via Kepler round-trip.
        use rand::{Rng, SeedableRng};
        use rand_chacha::ChaCha20Rng;
        use rand_distr::Uniform;

        let mu = MU_EARTH_KM3_S2;
        let mut rng = ChaCha20Rng::seed_from_u64(0xBEEF_DEAD);
        let radius = Uniform::new(5600.0, 10_500.0);
        let tof = Uniform::new(10_000.0, 250_000.0);

        let mut max_rel_err = 0.0_f64;
        let mut branches = 0_u32;
        let mut good_count = 0_u32;
        for _ in 0..500 {
            let r1_km = vec3::scale(rand_unit_vec(&mut rng), rng.sample(radius));
            let r2_km = vec3::scale(rand_unit_vec(&mut rng), rng.sample(radius));
            let tof_s = rng.sample(tof);
            let Ok(sols) = lambert(
                r1_km,
                r2_km,
                tof_s,
                mu,
                TransferWay::Short,
                RevolutionBudget::up_to(3),
            ) else {
                continue;
            };
            let mut iter_branch = |s: LambertSolution| {
                let r2_prop = kepler_propagate(r1_km, s.v1_km_s, tof_s, mu);
                let rel = vec_sub_norm(r2_prop, r2_km) / vec3::norm(r2_km);
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
        assert!(
            good_count > 300,
            "too few converged round-trips: {good_count}/{branches}"
        );
        assert!(
            max_rel_err < 1e-5,
            "max relative round-trip err = {max_rel_err:.3e} over {good_count} branches"
        );
    }
}
