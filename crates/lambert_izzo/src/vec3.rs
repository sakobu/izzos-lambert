//! Internal `[f64; 3]` vector helpers for the Lambert kernel.
//!
//! The crate's public surface uses plain `[f64; 3]` arrays so the core has
//! zero hard math-library dependencies. `nalgebra::Vector3<f64>` and
//! `glam::DVec3` already convert to/from `[f64; 3]` natively for callers
//! who want them.

// `num_traits::Float` provides `sqrt`/`sin`/`cos`/etc. on `f64` in `no_std`
// builds (via `libm`). In `std` builds the inherent methods take precedence
// and the trait import is dead but harmless — same pattern in every kernel
// module that touches transcendental math.
#[allow(unused_imports)]
use num_traits::Float as _;

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
pub(crate) fn norm(a: Vec3) -> f64 {
    dot(a, a).sqrt()
}

#[inline]
pub(crate) fn scale(a: Vec3, s: f64) -> Vec3 {
    [a[0] * s, a[1] * s, a[2] * s]
}
