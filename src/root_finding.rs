//! Householder + Halley iteration on `T(x, λ, M) = T_target`.
//!
//! Drives the search from the analytic initial guesses (Izzo Eqs. 30, 31)
//! to converged `(x, y)` pairs for every reachable revolution count.

use core::f64::consts::PI;

use crate::constants::{
    HALLEY_MAX_ITERS, HALLEY_TOL, HOUSEHOLDER_MAX_ITERS, HOUSEHOLDER_TOL_MULTI,
    HOUSEHOLDER_TOL_SINGLE,
};
use crate::error::LambertError;
use crate::geometry::Geometry;
use crate::tof::{compute_y, tof_derivatives, x_to_tof};

/// One converged root: `(x, y, n_revs, iters)`.
pub(crate) type XYTuple = (f64, f64, u32, u32);

/// Find every reachable `(x, y)` root for the given geometry and revolution
/// budget.
///
/// Returns single-rev first (always exists), then `(long-period, short-period)`
/// pairs for each `M` from `1` up to `min(multi_revs, ⌊T/π⌋)`.
///
/// # Errors
///
/// Propagates [`LambertError::NoConvergence`] / [`LambertError::SingularDenominator`]
/// from the underlying Householder iteration.
#[allow(clippy::similar_names)] // xl/xr, yl/yr, il/ir, x0l/x0r are long/short-period multi-rev branch pairs — Izzo §3.
pub(crate) fn find_xy(
    geom: &Geometry,
    multi_revs: u32,
) -> Result<Vec<XYTuple>, LambertError> {
    let lambda = geom.lambda;
    let big_t = geom.big_t;
    debug_assert!(lambda.abs() < 1.0);
    debug_assert!(big_t > 0.0);

    let t00 = lambda.acos() + lambda * (1.0 - lambda * lambda).sqrt();

    // Single revolution — always exists.
    let t1 = (2.0 / 3.0) * (1.0 - lambda * lambda * lambda);
    let x0 = initial_guess_single_rev(big_t, t00, t1, lambda);
    let (x, iters) = householder(x0, big_t, lambda, 0)?;
    let y = compute_y(x, lambda);
    let mut out = vec![(x, y, 0_u32, iters)];

    // Multi-rev branches. Iterate upward from M = 1; stop at the first
    // infeasible branch. T_min(M) is monotonically increasing in M, so once
    // a branch fails, all higher branches also fail. The Halley T_min check
    // only fires on the boundary branch where big_t lies below the analytic
    // minimum T(x=0, M) = t00 + M·π — at most one Halley call per call to
    // `find_xy`.
    for m in 1..=multi_revs {
        let m_pi = f64::from(m) * PI;

        // Quick reject: T_min(M) ≥ M·π always, so big_t < M·π ⇒ infeasible.
        if big_t < m_pi {
            break;
        }

        // Boundary regime: big_t below the analytic minimum at x = 0.
        // Confirm numerical T_min(M) ≤ big_t; otherwise drop the branch
        // (and all higher M, which have higher T_min).
        if big_t < t00 + m_pi {
            let (_x_min, t_min) = halley_t_min(lambda, m);
            if t_min > big_t {
                break;
            }
        }

        let (x0l, x0r) = initial_guess_multi_rev(big_t, m);
        let (xl, il) = householder(x0l, big_t, lambda, m)?;
        let yl = compute_y(xl, lambda);
        let (xr, ir) = householder(x0r, big_t, lambda, m)?;
        let yr = compute_y(xr, lambda);
        out.push((xl, yl, m, il));
        out.push((xr, yr, m, ir));
    }
    Ok(out)
}

/// Single-rev initial guess (Izzo Eq. 30 + the derivative-matched
/// hyperbolic starter just after).
fn initial_guess_single_rev(big_t: f64, t0: f64, t1: f64, lambda: f64) -> f64 {
    if big_t >= t0 {
        // T ≥ T0: high-energy / slow-transfer side.
        (t0 / big_t).powf(2.0 / 3.0) - 1.0
    } else if big_t <= t1 {
        // T ≤ T1: hyperbolic side; derivative-matched starter.
        2.5 * t1 * (t1 - big_t) / (big_t * (1.0 - lambda.powi(5))) + 1.0
    } else {
        // T1 < T < T0: log-linear interpolation in (ξ, τ) plane.
        // Paper Eq. 30 displays the middle branch as `(T0/T)^(log2(T1/T0)) − 1`,
        // which fails the boundary x = 1 at T = T1 (a known typo). The textual
        // derivation just before Eq. 30 (page 12) gives the correct form below;
        // PyKEP uses the same formulation.
        (t0 / big_t).powf(2_f64.ln() / (t0 / t1).ln()) - 1.0
    }
}

/// Multi-rev initial guesses (Izzo Eq. 31). Returns `(left, right)` —
/// long-period and short-period asymptotes respectively.
#[allow(clippy::similar_names)] // kl/kr, x0l/x0r are the long/short-period asymptote starters — Izzo Eq. 31.
fn initial_guess_multi_rev(big_t: f64, m: u32) -> (f64, f64) {
    let m_pi = f64::from(m) * PI;
    let kl = ((m_pi + PI) / (8.0 * big_t)).powf(2.0 / 3.0);
    let x0l = (kl - 1.0) / (kl + 1.0);
    let kr = ((8.0 * big_t) / m_pi).powf(2.0 / 3.0);
    let x0r = (kr - 1.0) / (kr + 1.0);
    (x0l, x0r)
}

/// Householder's 3rd-order iteration for `T(x) − T_target = 0`.
#[allow(clippy::similar_names)] // tol vs tof, and dt/ddt/dddt are 1st/2nd/3rd T-derivatives — Izzo Eq. 22.
fn householder(
    mut x: f64,
    big_t: f64,
    lambda: f64,
    m: u32,
) -> Result<(f64, u32), LambertError> {
    let tol = if m == 0 {
        HOUSEHOLDER_TOL_SINGLE
    } else {
        HOUSEHOLDER_TOL_MULTI
    };
    let mut last_step = f64::INFINITY;
    for it in 1..=HOUSEHOLDER_MAX_ITERS {
        let tof = x_to_tof(x, lambda, m);
        let (dt, ddt, dddt) = tof_derivatives(x, lambda, tof);
        let f = tof - big_t;
        let denom = dt * (dt * dt - f * ddt) + dddt * f * f / 6.0;
        if denom == 0.0 {
            return Err(LambertError::SingularDenominator { n_revs: m });
        }
        let delta = f * (dt * dt - 0.5 * f * ddt) / denom;
        x -= delta;
        last_step = delta.abs();
        if last_step < tol {
            return Ok((x, it));
        }
    }
    Err(LambertError::NoConvergence {
        iterations: HOUSEHOLDER_MAX_ITERS,
        last_step,
        n_revs: m,
    })
}

/// Locate the minimum of `T(x)` along the M-revolution branch starting at
/// `x = 0`. Halley iteration on `dT/dx = 0`. Returns `(x_min, T_min)`.
#[allow(clippy::similar_names)] // dt/ddt/dddt are 1st/2nd/3rd T-derivatives — Izzo Eq. 22.
fn halley_t_min(lambda: f64, m: u32) -> (f64, f64) {
    let mut x = 0.0;
    for _ in 0..HALLEY_MAX_ITERS {
        let tof = x_to_tof(x, lambda, m);
        let (dt, ddt, dddt) = tof_derivatives(x, lambda, tof);
        let denom = 2.0 * ddt * ddt - dt * dddt;
        if denom == 0.0 {
            break;
        }
        let delta = 2.0 * dt * ddt / denom;
        x -= delta;
        if delta.abs() < HALLEY_TOL {
            break;
        }
    }
    let t_min = x_to_tof(x, lambda, m);
    (x, t_min)
}
