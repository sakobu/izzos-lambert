//! Inline integration tests for `lambert_izzo`.
//!
//! Split per scenario for navigability; all submodules see the same crate
//! root via `crate::*`. Tests use realistic Earth/Sun SI values from
//! `lambert_izzo_test_support::bodies` — they are scenario-named, not
//! generic-math-named.

#![allow(clippy::similar_names, clippy::unwrap_used)] // r/r1/r2 test scenario inputs follow paper convention; tests are exempt from the lib's no-unwrap rule.

mod errors;
mod interop;
mod kepler_roundtrip;
mod multi_rev;
mod regimes;
mod single_rev;

use glam::DVec3;

/// Tolerance comparison for finite-precision scalar checks.
pub(crate) fn approx(a: f64, b: f64, tol: f64) -> bool {
    (a - b).abs() < tol
}

/// Euclidean distance between two `[f64; 3]` vectors via `glam`.
pub(crate) fn vec_sub_norm(a: [f64; 3], b: [f64; 3]) -> f64 {
    (DVec3::from_array(a) - DVec3::from_array(b)).length()
}
