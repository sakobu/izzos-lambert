//! Reference Kepler propagator used only for round-trip validation tests.
//!
//! Not part of the public API. Universal-variable formulation, solved by
//! Newton iteration on the universal Kepler equation.
//!
//! Inputs follow the same SI convention as the public Lambert API
//! (`r0_km`, `v0_km_s`, `dt_s`, `mu_km3_s2`).

// See `vec3.rs` for the rationale on this `allow(unused_imports)`.
#[allow(unused_imports)]
use num_traits::Float as _;

use crate::vec3::{self, Vec3};

/// Newton-iteration cap on the universal Kepler equation. Well above the
/// `<30` typical convergence count; pathological geometries hit it and bail
/// to NaN, which statistical tests filter explicitly.
const KEPLER_MAX_ITERS: u32 = 200;

/// Convergence tolerance on the universal-variable step `|Δχ|`.
/// Tight enough that the round-trip Kepler check resolves the Lambert
/// solver's `1e-9` / `1e-11` energy/momentum residuals.
const KEPLER_NEWTON_TOL: f64 = 1e-12;

/// Threshold below which the Newton denominator `r` is treated as zero
/// (would otherwise blow up `Δχ = -f/r`). Same magnitude as the convergence
/// tol — both are "essentially f64 noise on a unit-scale quantity".
const KEPLER_DENOM_EPS: f64 = 1e-12;

/// Divergence sentinel on `|χ|`. Universal-variable iteration on a poorly-
/// chosen initial guess can blow up; cap the magnitude to bail to NaN
/// instead of looping forever inside the iteration limit.
const KEPLER_CHI_DIVERGENCE: f64 = 1e12;

/// Branch threshold for the Stumpff `c2(ψ)`, `c3(ψ)` evaluator. Above this
/// magnitude, use the closed-form trig/hyperbolic expressions (numerically
/// fine far from `ψ = 0`); below it, use the Taylor series (closed forms
/// suffer catastrophic cancellation as `ψ → 0`).
const STUMPFF_SERIES_THRESHOLD: f64 = 1e-6;

/// Propagate `(r0_km, v0_km_s)` for `dt_s` under two-body gravity.
///
/// Returns a NaN-vector if the universal-variable Newton iteration diverges
/// (e.g. on near-parabolic geometries where the simple `χ₀ = √μ·dt·|α|`
/// initial guess is too poor to recover). Callers in statistical tests
/// filter divergent trials by checking componentwise finiteness.
#[must_use]
pub fn kepler_propagate(
    r0_km: [f64; 3],
    v0_km_s: [f64; 3],
    dt_s: f64,
    mu_km3_s2: f64,
) -> [f64; 3] {
    let r0n = vec3::norm(r0_km);
    let v0n2 = vec3::norm_squared(v0_km_s);
    let alpha = 2.0 / r0n - v0n2 / mu_km3_s2; // 1/a
    let sqrt_mu = mu_km3_s2.sqrt();
    let sigma0 = vec3::dot(r0_km, v0_km_s) / sqrt_mu;
    let mut chi = sqrt_mu * dt_s * alpha.abs(); // initial guess
    let nan_vec: Vec3 = [f64::NAN, f64::NAN, f64::NAN];
    let mut converged = false;
    for _ in 0..KEPLER_MAX_ITERS {
        let psi = alpha * chi * chi;
        let (c2, c3) = stumpff(psi);
        let r = chi * chi * c2 + sigma0 * chi * (1.0 - psi * c3) + r0n * (1.0 - psi * c2);
        if r.abs() < KEPLER_DENOM_EPS || !r.is_finite() {
            return nan_vec;
        }
        let f = sigma0 * chi * chi * c2 + (1.0 - alpha * r0n) * chi * chi * chi * c3 + r0n * chi
            - sqrt_mu * dt_s;
        let dchi = -f / r;
        chi += dchi;
        if !chi.is_finite() || chi.abs() > KEPLER_CHI_DIVERGENCE {
            return nan_vec;
        }
        if dchi.abs() < KEPLER_NEWTON_TOL {
            converged = true;
            break;
        }
    }
    if !converged {
        return nan_vec;
    }
    let psi = alpha * chi * chi;
    let (c2, c3) = stumpff(psi);
    let f = 1.0 - chi * chi / r0n * c2;
    let g = dt_s - chi * chi * chi / sqrt_mu * c3;
    vec3::add(vec3::scale(r0_km, f), vec3::scale(v0_km_s, g))
}

/// Stumpff functions `c2(ψ)`, `c3(ψ)` with a small-ψ series expansion.
fn stumpff(psi: f64) -> (f64, f64) {
    if psi > STUMPFF_SERIES_THRESHOLD {
        let s = psi.sqrt();
        ((1.0 - s.cos()) / psi, (s - s.sin()) / s.powi(3))
    } else if psi < -STUMPFF_SERIES_THRESHOLD {
        let s = (-psi).sqrt();
        ((1.0 - s.cosh()) / psi, (s.sinh() - s) / s.powi(3))
    } else {
        // Taylor coefficients for c2(ψ) = (1 − cos√ψ)/ψ and c3(ψ) = (√ψ − sin√ψ)/√ψ³
        // expanded about ψ = 0; denominators are factorials (2!, 4!, 6!) and
        // (3!, 5!, 7!) respectively.
        (
            0.5 - psi / 24.0 + psi * psi / 720.0,
            1.0 / 6.0 - psi / 120.0 + psi * psi / 5040.0,
        )
    }
}
