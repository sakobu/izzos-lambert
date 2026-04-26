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

use crate::constants::{
    BATTIN_THRESHOLD, HYPERGEOMETRIC_2F1_MAX_TERMS, HYPERGEOMETRIC_2F1_TOL, LAGRANGE_THRESHOLD,
};

/// Evaluate `T(x, λ, M)` using the regime-appropriate formulation.
pub(crate) fn x_to_tof(x: f64, lambda: f64, m: u32) -> f64 {
    let dist = (x - 1.0).abs();
    if dist <= BATTIN_THRESHOLD {
        x_to_tof_battin(x, lambda)
    } else if dist <= LAGRANGE_THRESHOLD {
        x_to_tof_lancaster(x, lambda, m)
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
        0.5 * a
            * a.sqrt()
            * ((alpha - alpha.sin()) - (beta - beta.sin()) + 2.0 * PI * f64::from(m))
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
fn x_to_tof_lancaster(x: f64, lambda: f64, m: u32) -> f64 {
    let y = compute_y(x, lambda);
    let one_minus_x2 = 1.0 - x * x;
    let psi = compute_psi(x, y, lambda);
    (((psi + f64::from(m) * PI) / one_minus_x2.abs().sqrt()) - x + lambda * y) / one_minus_x2
}

/// Auxiliary angle ψ (Izzo Eq. 17). `atan2` argument order is intentional —
/// Eq. 17 of the paper defines ψ as `atan2(numerator, denominator)` with
/// these specific terms.
fn compute_psi(x: f64, y: f64, lambda: f64) -> f64 {
    if x.abs() < 1.0 {
        // Elliptic regime, ψ ∈ [0, π].
        ((y - x * lambda) * (1.0 - x * x).sqrt()).atan2(x * y + lambda * (1.0 - x * x))
    } else {
        // Hyperbolic regime.
        let arg = (y - x * lambda) * (x * x - 1.0).sqrt();
        arg.asinh()
    }
}

/// Battin hypergeometric formulation (Izzo Eq. 20). Used near the parabolic
/// point `x = 1`. Implicitly `M = 0` — the near-parabolic regime is
/// single-revolution only.
fn x_to_tof_battin(x: f64, lambda: f64) -> f64 {
    let y = compute_y(x, lambda);
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
/// Consumed by the Householder iteration in [`crate::root_finding`].
#[allow(clippy::similar_names)] // dt/ddt/dddt are 1st/2nd/3rd derivatives; l2/l3/l5 are λ²/λ³/λ⁵ — Izzo Eq. 22.
pub(crate) fn tof_derivatives(x: f64, lambda: f64, tof: f64) -> (f64, f64, f64) {
    let y = compute_y(x, lambda);
    let one_m_x2 = 1.0 - x * x;
    let l2 = lambda * lambda;
    let l3 = l2 * lambda;
    let l5 = l3 * l2;

    let dt = (3.0 * tof * x - 2.0 + 2.0 * l3 * x / y) / one_m_x2;
    let ddt = (3.0 * tof + 5.0 * x * dt + 2.0 * (1.0 - l2) * l3 / (y * y * y)) / one_m_x2;
    let dddt =
        (7.0 * x * ddt + 8.0 * dt - 6.0 * (1.0 - l2) * l5 * x / (y.powi(5))) / one_m_x2;
    (dt, ddt, dddt)
}
