//! Internal `[f64; 3]` vector helpers for the Lambert kernel.
//!
//! The crate's public surface uses plain `[f64; 3]` arrays so the core has
//! zero hard math-library dependencies. Optional `nalgebra` and `glam`
//! features add `From`/`Into` conversions for ergonomic interop, but the
//! kernel itself stays in arrays.

pub(crate) type Vec3 = [f64; 3];

#[inline]
pub(crate) fn dot(a: Vec3, b: Vec3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
pub(crate) fn cross(a: Vec3, b: Vec3) -> Vec3 {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[inline]
pub(crate) fn norm_squared(a: Vec3) -> f64 {
    dot(a, a)
}

#[inline]
pub(crate) fn norm(a: Vec3) -> f64 {
    norm_squared(a).sqrt()
}

#[inline]
pub(crate) fn add(a: Vec3, b: Vec3) -> Vec3 {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

#[inline]
pub(crate) fn sub(a: Vec3, b: Vec3) -> Vec3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
pub(crate) fn scale(a: Vec3, s: f64) -> Vec3 {
    [a[0] * s, a[1] * s, a[2] * s]
}

#[inline]
pub(crate) fn normalize(a: Vec3) -> Vec3 {
    scale(a, 1.0 / norm(a))
}

/// Returns `Some(unit)` if `|a| >= threshold`, else `None`.
///
/// Mirrors `nalgebra::Vector3::try_normalize` so call sites read identically
/// to the previous implementation.
#[inline]
pub(crate) fn try_normalize(a: Vec3, threshold: f64) -> Option<Vec3> {
    let n = norm(a);
    if n >= threshold {
        Some(scale(a, 1.0 / n))
    } else {
        None
    }
}
