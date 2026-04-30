//! Time-of-flight evaluation along the Izzo `x → T(x, λ, M)` curve.
//!
//! Three regimes blend smoothly:
//! - Battin hypergeometric series (Izzo Eq. 20) for `|x − 1| ≤ BATTIN_THRESHOLD`.
//! - Lancaster–Blanchard (Eq. 18) for `BATTIN_THRESHOLD < |x − 1| ≤ LAGRANGE_THRESHOLD`.
//! - Lagrange (Eq. 9) elsewhere.
//!
//! Plus the analytic derivatives `dT/dx`, `d²T/dx²`, `d³T/dx³` (Eq. 22) used
//! by the Householder root finder in [`crate::root_finding`].

use core::f64::consts::PI;

// See `vec3.rs` for the rationale on this `allow(unused_imports)`.
#[allow(unused_imports)]
use num_traits::Float as _;

use crate::constants::{
    BATTIN_THRESHOLD, HYPERGEOMETRIC_2F1_MAX_TERMS, HYPERGEOMETRIC_2F1_TOL, LAGRANGE_THRESHOLD,
};

/// Evaluate `T(x, λ, M)` using the regime-appropriate formulation.
///
/// `y = sqrt(1 − λ²(1 − x²))` is supplied by the caller (already needed
/// by [`tof_derivatives_with_y`] on the next line of the Householder
/// loop) so we don't recompute it.
pub(crate) fn x_to_tof_with_y(x: f64, y: f64, lambda: f64, m: u32) -> f64 {
    let dist = (x - 1.0).abs();
    if m == 0 && dist <= BATTIN_THRESHOLD {
        x_to_tof_battin(x, y, lambda)
    } else if dist <= LAGRANGE_THRESHOLD {
        x_to_tof_lancaster(x, y, lambda, m)
    } else {
        x_to_tof_lagrange(x, lambda, m)
    }
}

/// `y = sqrt(1 − λ²(1 − x²))` (Izzo Eq. 7).
///
/// The `.max(0.0)` guards against `f64` round-off producing a tiny
/// negative argument when `|λ| → 1` and `x` is just past `1`. Same
/// idiom used in [`crate::geometry::Geometry::from_inputs`] for `λ` itself.
#[inline]
pub(crate) fn compute_y(x: f64, lambda: f64) -> f64 {
    (1.0 - lambda * lambda * (1.0 - x * x)).max(0.0).sqrt()
}

/// Lagrange formulation (Izzo Eq. 9). Stable away from the parabolic
/// point `x = 1`.
fn x_to_tof_lagrange(x: f64, lambda: f64, m: u32) -> f64 {
    let a = 1.0 / (1.0 - x * x);
    if a > 0.0 {
        // Ellipse.
        let alpha = 2.0 * x.acos();
        let mut beta = 2.0 * (lambda * lambda / a).abs().sqrt().asin();
        if lambda < 0.0 {
            beta = -beta;
        }
        0.5 * a * a.sqrt() * ((alpha - alpha.sin()) - (beta - beta.sin()) + 2.0 * PI * f64::from(m))
    } else {
        // Hyperbola — no multi-rev. From Eq. 15: λ·sqrt(x²−1) = sinh(β/2),
        // and (x²−1) = −1/a, so sinh(β/2) = λ·sqrt(−1/a).
        let alpha = 2.0 * x.acosh();
        let mut beta = 2.0 * (-lambda * lambda / a).sqrt().asinh();
        if lambda < 0.0 {
            beta = -beta;
        }
        -0.5 * a * (-a).sqrt() * ((beta - beta.sinh()) - (alpha - alpha.sinh()))
    }
}

/// Lancaster–Blanchard formulation (Izzo Eq. 18). Stable in the band
/// `BATTIN_THRESHOLD < |x − 1| ≤ LAGRANGE_THRESHOLD`.
fn x_to_tof_lancaster(x: f64, y: f64, lambda: f64, m: u32) -> f64 {
    let one_m_x2 = 1.0 - x * x;
    let psi = compute_psi(x, y, lambda, one_m_x2);
    (((psi + f64::from(m) * PI) / one_m_x2.abs().sqrt()) - x + lambda * y) / one_m_x2
}

/// Auxiliary angle ψ (Izzo Eq. 17). `atan2` argument order is intentional —
/// Eq. 17 of the paper defines ψ as `atan2(numerator, denominator)` with
/// these specific terms. `one_m_x2 = 1 − x²` is supplied by the caller
/// (already needed for the Lancaster denominator) so we don't recompute it.
fn compute_psi(x: f64, y: f64, lambda: f64, one_m_x2: f64) -> f64 {
    if x.abs() < 1.0 {
        // Elliptic regime, ψ ∈ [0, π].
        ((y - x * lambda) * one_m_x2.sqrt()).atan2(x * y + lambda * one_m_x2)
    } else {
        // Hyperbolic regime.
        ((y - x * lambda) * (-one_m_x2).sqrt()).asinh()
    }
}

/// Battin hypergeometric formulation (Izzo Eq. 20). Used near the parabolic
/// point `x = 1`. Implicitly `M = 0` — the near-parabolic regime is
/// single-revolution only.
fn x_to_tof_battin(x: f64, y: f64, lambda: f64) -> f64 {
    let eta = y - lambda * x;
    let s1 = 0.5 * (1.0 - lambda - x * eta);
    let q = (4.0 / 3.0) * hypergeometric_2f1_special(s1);
    0.5 * (eta * eta * eta * q + 4.0 * lambda * eta)
}

/// Direct series for `2F1(3, 1; 5/2; z)`. Converges fast for `|z| < 1`.
///
/// On hitting [`HYPERGEOMETRIC_2F1_MAX_TERMS`], returns the partial sum —
/// the cap is a safety guard against pathological inputs (NaN), not an
/// expected exit. In normal use the series converges in `< 30` terms.
fn hypergeometric_2f1_special(z: f64) -> f64 {
    let mut s = 1.0;
    let mut c = 1.0;
    for j in 0..HYPERGEOMETRIC_2F1_MAX_TERMS {
        let jf = f64::from(j);
        c *= (3.0 + jf) * (1.0 + jf) / (2.5 + jf) * z / (jf + 1.0);
        let s_new = s + c;
        if (s_new - s).abs() < HYPERGEOMETRIC_2F1_TOL {
            return s_new;
        }
        s = s_new;
    }
    s
}

/// Analytic derivatives `dT/dx`, `d²T/dx²`, `d³T/dx³` (Izzo Eq. 22).
/// Consumed by the Householder iteration in [`crate::root_finding`]. `y`
/// is supplied by the caller — see [`x_to_tof_with_y`] for the rationale.
#[allow(clippy::similar_names)] // dt/ddt/dddt are 1st/2nd/3rd derivatives; l2/l3/l5 are λ²/λ³/λ⁵ — Izzo Eq. 22.
pub(crate) fn tof_derivatives_with_y(
    x: f64,
    y: f64,
    lambda: f64,
    tof: f64,
) -> (f64, f64, f64) {
    let one_m_x2 = 1.0 - x * x;
    let l2 = lambda * lambda;
    let l3 = l2 * lambda;
    let l5 = l3 * l2;

    let dt = (3.0 * tof * x - 2.0 + 2.0 * l3 * x / y) / one_m_x2;
    let ddt = (3.0 * tof + 5.0 * x * dt + 2.0 * (1.0 - l2) * l3 / (y * y * y)) / one_m_x2;
    let dddt = (7.0 * x * ddt + 8.0 * dt - 6.0 * (1.0 - l2) * l5 * x / (y.powi(5))) / one_m_x2;
    (dt, ddt, dddt)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn multi_rev_near_parabolic_keeps_revolution_term() {
        let x = 0.995;
        let lambda = 0.4;
        let m = 1;

        let y = compute_y(x, lambda);
        let got = x_to_tof_with_y(x, y, lambda, m);
        let expected = x_to_tof_lancaster(x, y, lambda, m);

        assert!(
            (got - expected).abs() < 1e-12,
            "multi-rev TOF dispatch dropped M term: got {got}, expected {expected}"
        );
    }

    /// Cross-check the analytic Eq. 22 derivatives against centered finite
    /// differences. Catches algebra slips immediately rather than waiting
    /// for the Householder convergence stats to regress.
    #[test]
    #[allow(clippy::similar_names)] // t_p1/t_m1/t_p2/t_m2 are forward/backward TOF samples at ±h, ±2h;
    // dt_fd/ddt_fd/dddt_fd are FD-derived 1st/2nd/3rd derivatives matching tof_derivatives_with_y's tuple — Izzo Eq. 22.
    fn analytic_derivatives_match_finite_differences() {
        // Sample points chosen to span the Lagrange (large `|x − 1|`) and
        // Lancaster (mid-band) regimes, single- and multi-rev.
        // The Battin regime (`|x − 1| ≤ 0.01`) is excluded — its TOF is
        // computed from a different formula (Eq. 20) and the Householder
        // loop only feeds derivatives through Eq. 22 outside that band.
        let cases = [
            (0.5_f64, 0.3_f64, 0_u32),
            (-0.3, 0.5, 0),
            (0.85, 0.2, 0),
            (0.5, 0.3, 2),
            (-0.5, -0.4, 1),
            (0.3, 0.7, 3),
        ];

        let h = 1e-3;
        let tof_at = |x: f64, lambda: f64, m: u32| {
            x_to_tof_with_y(x, compute_y(x, lambda), lambda, m)
        };

        for (x, lambda, m) in cases {
            let t = tof_at(x, lambda, m);
            let y = compute_y(x, lambda);
            let (dt_a, ddt_a, dddt_a) = tof_derivatives_with_y(x, y, lambda, t);

            let t_p1 = tof_at(x + h, lambda, m);
            let t_m1 = tof_at(x - h, lambda, m);
            let t_p2 = tof_at(x + 2.0 * h, lambda, m);
            let t_m2 = tof_at(x - 2.0 * h, lambda, m);

            let dt_fd = (t_p1 - t_m1) / (2.0 * h);
            let ddt_fd = (t_p1 - 2.0 * t + t_m1) / (h * h);
            let dddt_fd = (t_p2 - 2.0 * t_p1 + 2.0 * t_m1 - t_m2) / (2.0 * h * h * h);

            let rel_err = |a: f64, b: f64| (a - b).abs() / a.abs().max(1.0);

            assert!(
                rel_err(dt_a, dt_fd) < 1e-5,
                "dT/dx mismatch at (x={x}, λ={lambda}, M={m}): \
                 analytic = {dt_a}, FD = {dt_fd}",
            );
            assert!(
                rel_err(ddt_a, ddt_fd) < 1e-3,
                "d²T/dx² mismatch at (x={x}, λ={lambda}, M={m}): \
                 analytic = {ddt_a}, FD = {ddt_fd}",
            );
            assert!(
                rel_err(dddt_a, dddt_fd) < 1e-1,
                "d³T/dx³ mismatch at (x={x}, λ={lambda}, M={m}): \
                 analytic = {dddt_a}, FD = {dddt_fd}",
            );
        }
    }
}
