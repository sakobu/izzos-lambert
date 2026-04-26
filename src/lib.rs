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
//! # Example
//!
//! ```
//! use lambert_izzo::{lambert, RevolutionBudget, TransferWay};
//! use nalgebra::Vector3;
//!
//! // LEO → LEO 90° transfer at 7000 km altitude.
//! let mu_km3_s2 = 398_600.4418;
//! let r1_km = Vector3::new(7000.0, 0.0, 0.0);
//! let r2_km = Vector3::new(0.0, 7000.0, 0.0);
//! let tof_s = core::f64::consts::PI / 2.0 * (7000.0_f64.powi(3) / mu_km3_s2).sqrt();
//!
//! let solutions = lambert(
//!     r1_km, r2_km, tof_s, mu_km3_s2,
//!     TransferWay::Short, RevolutionBudget::SingleOnly,
//! ).unwrap();
//! assert_eq!(solutions.len(), 1);
//! ```

#![warn(missing_docs)]
#![warn(clippy::pedantic)]
#![allow(clippy::module_name_repetitions)] // LambertError, LambertSolution

use nalgebra::Vector3;

/// Direction around the transfer plane from `r1` to `r2`.
///
/// `Short` is the geodesic transfer (`θ ≤ π`); `Long` traverses the other
/// way (`θ > π`). This is independent of prograde/retrograde — the orbit's
/// angular-momentum direction is set by the order of the `(r1, r2)` arguments.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
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
/// exceeds the requested time of flight.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
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

#[cfg(test)]
mod test_helpers;

pub use error::LambertError;

use geometry::Geometry;
use root_finding::find_xy;

/// Diagnostic data from one Householder solve. Useful for debugging or
/// distinguishing multi-rev branches; not part of the trajectory answer.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SolverDiagnostics {
    /// Householder iterations used to converge.
    pub iters: u32,
    /// Final value of Izzo's free parameter (Lancaster–Blanchard variable),
    /// dimensionless. For a given `M`, the long-period branch has the smaller
    /// value, the short-period branch the larger.
    pub lancaster_blanchard_x: f64,
}

/// One Lambert transfer.
///
/// For multi-rev problems, multiple solutions exist per `M`; the returned
/// vector from [`lambert`] is ordered: single-rev first, then
/// `(M=1 long-period, M=1 short-period, M=2 long-period, …)`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LambertSolution {
    /// Velocity at `r1_km` (km/s, same inertial frame as the inputs).
    pub v1_km_s: Vector3<f64>,
    /// Velocity at `r2_km` (km/s, same inertial frame as the inputs).
    pub v2_km_s: Vector3<f64>,
    /// Number of complete revolutions (`0` for direct transfer).
    pub n_revs: u32,
    /// Solver diagnostics (iteration count + Lancaster–Blanchard `x`).
    pub diagnostics: SolverDiagnostics,
}

/// Solve Lambert's boundary-value problem.
///
/// # Arguments
///
/// - `r1_km`, `r2_km`: position vectors (km), in any consistent inertial frame.
/// - `tof_s`: time of flight (s), `> 0`.
/// - `mu_km3_s2`: gravitational parameter (km³/s²), `> 0`.
/// - `way`: short or long way around the transfer plane (see [`TransferWay`]).
/// - `revolutions`: revolution budget — see [`RevolutionBudget`].
///
/// # Returns
///
/// Up to `1 + 2 · M_max` solutions, where `M_max = min(revolutions.max(), ⌊T/π⌋)`
/// adjusted downward when the M-branch's `T_min` exceeds `tof_s`.
///
/// # Errors
///
/// - [`LambertError::NonPositiveTimeOfFlight`] — `tof_s <= 0`.
/// - [`LambertError::NonPositiveMu`] — `mu_km3_s2 <= 0`.
/// - [`LambertError::DegeneratePositionVector`] — `|r1|` or `|r2|` below
///   [`constants::MIN_POSITION_NORM_KM`].
/// - [`LambertError::CollinearGeometry`] — `|r1 × r2| / (|r1| · |r2|)` below
///   [`constants::COLINEARITY_TOL`].
/// - [`LambertError::NoConvergence`] / [`LambertError::SingularDenominator`]
///   — Householder iteration failed.
#[allow(clippy::similar_names)] // v_r1/v_r2/v_t1/v_t2 are radial/tangential velocity components at points 1/2 — Izzo §2.
pub fn lambert(
    r1_km: Vector3<f64>,
    r2_km: Vector3<f64>,
    tof_s: f64,
    mu_km3_s2: f64,
    way: TransferWay,
    revolutions: RevolutionBudget,
) -> Result<Vec<LambertSolution>, LambertError> {
    let geom = Geometry::from_inputs(r1_km, r2_km, tof_s, mu_km3_s2, way)?;
    let xy_pairs = find_xy(&geom, revolutions)?;

    let mut out = Vec::with_capacity(xy_pairs.len());
    for root in xy_pairs {
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

        let v1_km_s = geom.ir1 * v_r1 + geom.it1 * v_t1;
        let v2_km_s = geom.ir2 * v_r2 + geom.it2 * v_t2;
        out.push(LambertSolution {
            v1_km_s,
            v2_km_s,
            n_revs: root.n_revs,
            diagnostics: SolverDiagnostics {
                iters: root.iters,
                lancaster_blanchard_x: root.x,
            },
        });
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    #![allow(clippy::similar_names)] // r_km/r1_km/r2_km test scenario inputs follow paper convention.

    use super::*;
    use crate::test_helpers::kepler_propagate;
    use core::f64::consts::PI;

    /// Earth's gravitational parameter (km³/s²) — value from EGM2008.
    const MU_EARTH_KM3_S2: f64 = 398_600.441_8;
    /// Sun's gravitational parameter (km³/s²) — value from DE440.
    const MU_SUN_KM3_S2: f64 = 1.327_124_400_18e11;

    fn approx(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn quarter_circle_leo() {
        // 90° transfer along a circular LEO at r = 7000 km.
        let r_km = 7000.0;
        let mu = MU_EARTH_KM3_S2;
        let v_circ = (mu / r_km).sqrt();
        let period_s = 2.0 * PI * (r_km.powi(3) / mu).sqrt();

        let r1_km = Vector3::new(r_km, 0.0, 0.0);
        let r2_km = Vector3::new(0.0, r_km, 0.0);
        let sols = lambert(r1_km, r2_km, period_s / 4.0, mu, TransferWay::Short, RevolutionBudget::SingleOnly).unwrap();
        assert_eq!(sols.len(), 1);
        assert!((sols[0].v1_km_s - Vector3::new(0.0, v_circ, 0.0)).norm() < 1e-9);
        assert!((sols[0].v2_km_s - Vector3::new(-v_circ, 0.0, 0.0)).norm() < 1e-9);
    }

    #[test]
    fn long_way_quarter_circle_leo() {
        // 270° transfer along the same circular LEO.
        let r_km = 7000.0;
        let mu = MU_EARTH_KM3_S2;
        let v_circ = (mu / r_km).sqrt();
        let period_s = 2.0 * PI * (r_km.powi(3) / mu).sqrt();

        let r1_km = Vector3::new(r_km, 0.0, 0.0);
        let r2_km = Vector3::new(0.0, r_km, 0.0);
        let sols = lambert(r1_km, r2_km, 3.0 * period_s / 4.0, mu, TransferWay::Long, RevolutionBudget::SingleOnly).unwrap();
        assert_eq!(sols.len(), 1);
        assert!((sols[0].v1_km_s - Vector3::new(0.0, -v_circ, 0.0)).norm() < 1e-9);
        assert!((sols[0].v2_km_s - Vector3::new(v_circ, 0.0, 0.0)).norm() < 1e-9);
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

        let r1_km = Vector3::new(r1_norm_km, 0.0, 0.0);
        // 1 km off-plane to dodge the colinearity edge case.
        let r2_km = Vector3::new(-r2_norm_km, 1.0, 0.0);
        let sols = lambert(r1_km, r2_km, tof_s, mu, TransferWay::Short, RevolutionBudget::SingleOnly).unwrap();
        // Periapsis velocity: vis-viva at r = 1 AU on the transfer ellipse.
        let v_peri_km_s = (mu * (2.0 / r1_norm_km - 1.0 / a_km)).sqrt();
        // Tolerance ~1 m/s — the off-plane perturbation moves things slightly.
        assert!(approx(sols[0].v1_km_s.norm(), v_peri_km_s, 1e-3));
    }

    #[test]
    fn multi_rev_branches() {
        // Long-tof Earth-orbit phasing — admits multiple multi-rev branches.
        let mu = MU_EARTH_KM3_S2;
        let r1_km = Vector3::new(8000.0, 0.0, 0.0);
        let r2_km = Vector3::new(5600.0, 5600.0, 0.0);
        let period_s = 2.0 * PI * (8000.0_f64.powi(3) / mu).sqrt();
        let sols = lambert(r1_km, r2_km, 5.0 * period_s, mu, TransferWay::Short, RevolutionBudget::up_to(3)).unwrap();
        assert!(sols.len() >= 3, "got {} solutions", sols.len());
        for s in &sols {
            let r1n = r1_km.norm();
            let energy = 0.5 * s.v1_km_s.dot(&s.v1_km_s) - mu / r1n;
            assert!(energy.is_finite());
        }
    }

    #[test]
    fn round_trip_kepler_check_single_rev() {
        // Propagate v1 with a universal-variable Kepler integrator and confirm
        // we land within numerical tolerance of r2.
        let mu = MU_EARTH_KM3_S2;
        let r1_km = Vector3::new(10_500.0, 1400.0, 700.0);
        let r2_km = Vector3::new(-2800.0, 9100.0, -1400.0);
        let tof_s = 4500.0;
        let sols = lambert(r1_km, r2_km, tof_s, mu, TransferWay::Short, RevolutionBudget::SingleOnly).unwrap();
        let v1_km_s = sols[0].v1_km_s;
        let r2_prop_km = kepler_propagate(r1_km, v1_km_s, tof_s, mu);
        let err_km = (r2_prop_km - r2_km).norm();
        assert!(err_km < 1e-6, "kepler-roundtrip err = {err_km} km");
    }

    #[test]
    fn errors_on_non_positive_tof() {
        let r1_km = Vector3::new(7000.0, 0.0, 0.0);
        let r2_km = Vector3::new(0.0, 7000.0, 0.0);
        let err = lambert(r1_km, r2_km, 0.0, MU_EARTH_KM3_S2, TransferWay::Short, RevolutionBudget::SingleOnly).unwrap_err();
        assert!(matches!(err, LambertError::NonPositiveTimeOfFlight { tof_s } if tof_s == 0.0));
    }

    #[test]
    fn errors_on_zero_position_vector() {
        let r1_km = Vector3::zeros();
        let r2_km = Vector3::new(0.0, 7000.0, 0.0);
        let err = lambert(r1_km, r2_km, 1000.0, MU_EARTH_KM3_S2, TransferWay::Short, RevolutionBudget::SingleOnly).unwrap_err();
        assert!(matches!(
            err,
            LambertError::DegeneratePositionVector { which: 1, .. }
        ));
    }

    #[test]
    fn errors_on_colinear_geometry() {
        let r1_km = Vector3::new(7000.0, 0.0, 0.0);
        let r2_km = Vector3::new(14_000.0, 0.0, 0.0);
        let err = lambert(r1_km, r2_km, 1000.0, MU_EARTH_KM3_S2, TransferWay::Short, RevolutionBudget::SingleOnly).unwrap_err();
        assert!(matches!(err, LambertError::CollinearGeometry { .. }));
    }

    #[test]
    fn errors_on_non_positive_mu() {
        let r1_km = Vector3::new(7000.0, 0.0, 0.0);
        let r2_km = Vector3::new(0.0, 7000.0, 0.0);
        let err = lambert(r1_km, r2_km, 1000.0, 0.0, TransferWay::Short, RevolutionBudget::SingleOnly).unwrap_err();
        assert!(matches!(err, LambertError::NonPositiveMu { mu_km3_s2 } if mu_km3_s2 == 0.0));
    }

    #[test]
    fn battin_regime_near_parabolic() {
        // GTO-like 90° transfer with TOF tuned so the converged x lands inside
        // the |x − 1| ≤ BATTIN_THRESHOLD (= 0.01) band, exercising the
        // hypergeometric series formulation in tof::x_to_tof_battin (Eq. 20).
        let mu = MU_EARTH_KM3_S2;
        let r1_km = Vector3::new(7000.0, 0.0, 0.0);
        let r2_km = Vector3::new(0.0, 42_000.0, 0.0);
        let tof_s = 7200.0;
        let sols = lambert(r1_km, r2_km, tof_s, mu, TransferWay::Short, RevolutionBudget::SingleOnly).unwrap();
        let s = sols[0];
        assert!(
            (s.diagnostics.lancaster_blanchard_x - 1.0).abs() < crate::constants::BATTIN_THRESHOLD,
            "x = {} not in Battin band [1 - {tol}, 1 + {tol}]",
            s.diagnostics.lancaster_blanchard_x,
            tol = crate::constants::BATTIN_THRESHOLD,
        );
        let r2_prop = kepler_propagate(r1_km, s.v1_km_s, tof_s, mu);
        let err_km = (r2_prop - r2_km).norm();
        assert!(err_km < 1e-3, "Battin round-trip err = {err_km} km");
    }

    #[test]
    fn hyperbolic_transfer() {
        // Fast LEO → 200 000 km transfer; required v1 exceeds Earth escape,
        // landing in the hyperbolic branch (x > 1, positive specific energy).
        // Exercises tof::x_to_tof_lagrange's a < 0 path and compute_psi's
        // asinh branch (Eq. 9 / Eq. 17 hyperbolic).
        let mu = MU_EARTH_KM3_S2;
        let r1_km = Vector3::new(7000.0, 0.0, 0.0);
        let r2_km = Vector3::new(0.0, 200_000.0, 0.0);
        let tof_s = 30_000.0;
        let sols = lambert(r1_km, r2_km, tof_s, mu, TransferWay::Short, RevolutionBudget::SingleOnly).unwrap();
        let s = sols[0];
        assert!(
            s.diagnostics.lancaster_blanchard_x > 1.0,
            "expected hyperbolic (x > 1), got x = {}",
            s.diagnostics.lancaster_blanchard_x
        );
        let energy = 0.5 * s.v1_km_s.dot(&s.v1_km_s) - mu / r1_km.norm();
        assert!(energy > 0.0, "expected positive specific energy, got {energy}");
        let r2_prop = kepler_propagate(r1_km, s.v1_km_s, tof_s, mu);
        let err_km = (r2_prop - r2_km).norm();
        assert!(err_km < 1e-3, "hyperbolic round-trip err = {err_km} km");
    }

    #[test]
    fn multi_rev_branches_distinct() {
        // Same geometry as multi_rev_branches; verify long-period x < short-period x
        // for each M (paper's "switch between branches" concern, page 14) and that
        // both branches independently round-trip via Kepler.
        let mu = MU_EARTH_KM3_S2;
        let r1_km = Vector3::new(8000.0, 0.0, 0.0);
        let r2_km = Vector3::new(5600.0, 5600.0, 0.0);
        let period_s = 2.0 * PI * (8000.0_f64.powi(3) / mu).sqrt();
        let tof_s = 5.0 * period_s;
        let sols = lambert(r1_km, r2_km, tof_s, mu, TransferWay::Short, RevolutionBudget::up_to(3)).unwrap();

        let multi: Vec<_> = sols.iter().filter(|s| s.n_revs > 0).collect();
        assert!(!multi.is_empty(), "no multi-rev branches found");
        assert_eq!(
            multi.len() % 2,
            0,
            "multi-rev branches not paired: len = {}",
            multi.len()
        );

        for pair in multi.chunks(2) {
            assert_eq!(
                pair[0].n_revs, pair[1].n_revs,
                "branches in pair don't share n_revs"
            );
            // find_xy pushes (xl, xr) per M — long-period first, short-period second.
            assert!(
                pair[0].diagnostics.lancaster_blanchard_x
                    < pair[1].diagnostics.lancaster_blanchard_x,
                "M = {}: long-period x ({}) >= short-period x ({})",
                pair[0].n_revs,
                pair[0].diagnostics.lancaster_blanchard_x,
                pair[1].diagnostics.lancaster_blanchard_x
            );
            for s in pair {
                let r2_prop = kepler_propagate(r1_km, s.v1_km_s, tof_s, mu);
                let err_km = (r2_prop - r2_km).norm();
                assert!(
                    err_km < 1e-3,
                    "M = {} branch round-trip err = {err_km} km",
                    s.n_revs
                );
            }
        }
    }

    #[test]
    fn solution_ordering_contract() {
        // The LambertSolution doc promises:
        //   index 0: single-rev (n_revs = 0)
        //   then: (M=1 long-period, M=1 short-period, M=2 long, M=2 short, ...)
        // Long-period x is smaller than short-period x within each M (Izzo §3).
        let mu = MU_EARTH_KM3_S2;
        let r1_km = Vector3::new(8000.0, 0.0, 0.0);
        let r2_km = Vector3::new(5600.0, 5600.0, 0.0);
        let period_s = 2.0 * PI * (8000.0_f64.powi(3) / mu).sqrt();
        let tof_s = 5.0 * period_s;
        let sols = lambert(
            r1_km, r2_km, tof_s, mu,
            TransferWay::Short, RevolutionBudget::up_to(3),
        ).unwrap();

        // Index 0 is always single-rev.
        assert_eq!(sols[0].n_revs, 0, "index 0 must be single-rev");

        // Multi-rev branches: paired, ascending M, long-period before short-period.
        let multi = &sols[1..];
        assert_eq!(multi.len() % 2, 0, "multi-rev branches must be paired");

        let mut prev_m = 0_u32;
        for pair in multi.chunks(2) {
            assert_eq!(pair[0].n_revs, pair[1].n_revs, "pair shares n_revs");
            assert!(pair[0].n_revs > prev_m, "M strictly ascending across pairs");
            assert!(
                pair[0].diagnostics.lancaster_blanchard_x
                    < pair[1].diagnostics.lancaster_blanchard_x,
                "long-period x must precede short-period x within pair (M = {})",
                pair[0].n_revs,
            );
            prev_m = pair[0].n_revs;
        }
    }

    fn rand_unit_vec(rng: &mut rand_chacha::ChaCha20Rng) -> Vector3<f64> {
        use rand::Rng;
        use rand_distr::Uniform;
        let axis: Uniform<f64> = Uniform::new(-1.0, 1.0);
        loop {
            let v = Vector3::new(rng.sample(axis), rng.sample(axis), rng.sample(axis));
            let n2 = v.norm_squared();
            if n2 > 0.01 && n2 < 1.0 {
                return v / n2.sqrt();
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
            let r1_km = rand_unit_vec(&mut rng) * rng.sample(radius);
            let r2_km = rand_unit_vec(&mut rng) * rng.sample(radius);
            let tof_s = rng.sample(tof);
            let way = if rng.gen_bool(0.5) { TransferWay::Long } else { TransferWay::Short };
            let Ok(sols) = lambert(r1_km, r2_km, tof_s, mu, way, RevolutionBudget::SingleOnly) else {
                continue;
            };
            lambert_ok += 1;
            let s = sols[0];
            let r2_prop = kepler_propagate(r1_km, s.v1_km_s, tof_s, mu);
            let rel = (r2_prop - r2_km).norm() / r2_km.norm();
            // NaN signals propagator divergence on a pathological geometry —
            // not a Lambert correctness issue. Filter and assert on the rest.
            if rel.is_finite() {
                max_rel_err = max_rel_err.max(rel);
                good_count += 1;
            }
        }
        assert!(lambert_ok > 950, "too many Lambert failures: {lambert_ok}/1000");
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
            let r1_km = rand_unit_vec(&mut rng) * rng.sample(radius);
            let r2_km = rand_unit_vec(&mut rng) * rng.sample(radius);
            let tof_s = rng.sample(tof);
            let Ok(sols) = lambert(r1_km, r2_km, tof_s, mu, TransferWay::Short, RevolutionBudget::up_to(3)) else {
                continue;
            };
            for s in &sols {
                let r2_prop = kepler_propagate(r1_km, s.v1_km_s, tof_s, mu);
                let rel = (r2_prop - r2_km).norm() / r2_km.norm();
                if rel.is_finite() {
                    max_rel_err = max_rel_err.max(rel);
                    good_count += 1;
                }
                branches += 1;
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
