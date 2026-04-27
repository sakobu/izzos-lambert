//! Lambert problem geometry: chord, semi-perimeter, λ, and transfer-plane basis.
//!
//! Constructed once per [`crate::lambert`] call and threaded through the
//! root-finding and velocity-reconstruction stages — eliminates re-derivation
//! drift across modules.

use crate::constants::{COLINEARITY_TOL, MIN_POSITION_NORM_KM};
use crate::error::LambertError;
use crate::vec3::{self, Vec3};

/// Pre-computed geometry for a Lambert boundary problem.
///
/// All scalars and unit vectors derive from `(r1_km, r2_km, tof_s, mu_km3_s2,
/// way)`. The solver kernels consume `lambda` and `big_t` (non-dimensional
/// TOF); the velocity reconstruction in [`crate::lambert`] consumes the rest.
///
/// Field names are math-domain names (paper symbols), not unit-tagged —
/// these are crate-private intermediates, not public state.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Geometry {
    /// Izzo's λ parameter, sign-corrected for short/long way (Eq. 7).
    pub lambda: f64,
    /// Non-dimensional time of flight `T = sqrt(2μ / s³) · tof`.
    pub big_t: f64,
    /// `sqrt(μ · s / 2)` — velocity scale for reconstruction.
    pub gamma: f64,
    /// `(r1n − r2n) / c`.
    pub rho: f64,
    /// `sqrt(1 − ρ²)`.
    pub sigma: f64,
    /// `|r1|`.
    pub r1n: f64,
    /// `|r2|`.
    pub r2n: f64,
    /// Unit vector along `r1`.
    pub ir1: Vec3,
    /// Unit vector along `r2`.
    pub ir2: Vec3,
    /// In-plane tangent at `r1` (sign-corrected for long-way transfers).
    pub it1: Vec3,
    /// In-plane tangent at `r2` (sign-corrected for long-way transfers).
    pub it2: Vec3,
}

impl Geometry {
    /// Build the geometry from raw inputs, validating scalars and the
    /// transfer plane along the way.
    ///
    /// # Errors
    ///
    /// - [`LambertError::NonFiniteInput`] — any public scalar input or position
    ///   vector component is `NaN`, `+inf`, or `-inf`.
    /// - [`LambertError::NonPositiveTimeOfFlight`] — `tof_s <= 0`.
    /// - [`LambertError::NonPositiveMu`] — `mu_km3_s2 <= 0`.
    /// - [`LambertError::DegeneratePositionVector`] — `|r1|` or `|r2|`
    ///   below [`MIN_POSITION_NORM_KM`].
    /// - [`LambertError::CollinearGeometry`] — `|r1 × r2| / (|r1| · |r2|)`
    ///   below [`COLINEARITY_TOL`].
    #[allow(clippy::similar_names)] // ir1/ir2 are radial unit vectors, it1/it2 tangential — Izzo Eq. 5–7.
    pub(crate) fn from_inputs(
        r1_km: Vec3,
        r2_km: Vec3,
        tof_s: f64,
        mu_km3_s2: f64,
        way: crate::TransferWay,
    ) -> Result<Self, LambertError> {
        validate_finite_vector("r1_km", r1_km)?;
        validate_finite_vector("r2_km", r2_km)?;
        validate_finite_scalar("tof_s", tof_s)?;
        validate_finite_scalar("mu_km3_s2", mu_km3_s2)?;

        if tof_s <= 0.0 {
            return Err(LambertError::NonPositiveTimeOfFlight { tof_s });
        }
        if mu_km3_s2 <= 0.0 {
            return Err(LambertError::NonPositiveMu { mu_km3_s2 });
        }

        let r1n = vec3::norm(r1_km);
        let r2n = vec3::norm(r2_km);
        if r1n < MIN_POSITION_NORM_KM {
            return Err(LambertError::DegeneratePositionVector {
                which: 1,
                norm_km: r1n,
            });
        }
        if r2n < MIN_POSITION_NORM_KM {
            return Err(LambertError::DegeneratePositionVector {
                which: 2,
                norm_km: r2n,
            });
        }

        let chord = vec3::sub(r2_km, r1_km);
        let c = vec3::norm(chord);
        let s = 0.5 * (r1n + r2n + c);

        let ir1 = vec3::scale(r1_km, 1.0 / r1n);
        let ir2 = vec3::scale(r2_km, 1.0 / r2n);
        let ih_raw = vec3::cross(ir1, ir2);
        let sin_angle = vec3::norm(ih_raw);
        let Some(ih) = vec3::try_normalize(ih_raw, COLINEARITY_TOL) else {
            return Err(LambertError::CollinearGeometry { sin_angle });
        };

        // λ² = 1 − c/s, λ ∈ [-1, 1]. Default convention: prograde
        // (counter-clockwise about ih) with θ ∈ [0, π] → λ > 0.
        // For long-way we flip λ AND the tangent direction.
        let mut lambda = (1.0 - c / s).max(0.0).sqrt();
        let (it1_raw, it2_raw) = if matches!(way, crate::TransferWay::Long) {
            lambda = -lambda;
            (vec3::cross(ir1, ih), vec3::cross(ir2, ih))
        } else {
            (vec3::cross(ih, ir1), vec3::cross(ih, ir2))
        };
        let it1 = vec3::normalize(it1_raw);
        let it2 = vec3::normalize(it2_raw);

        let big_t = (2.0 * mu_km3_s2 / (s * s * s)).sqrt() * tof_s;
        let gamma = (mu_km3_s2 * s / 2.0).sqrt();
        let rho = (r1n - r2n) / c;
        let sigma = (1.0 - rho * rho).max(0.0).sqrt();

        Ok(Self {
            lambda,
            big_t,
            gamma,
            rho,
            sigma,
            r1n,
            r2n,
            ir1,
            ir2,
            it1,
            it2,
        })
    }
}

fn validate_finite_scalar(parameter: &'static str, value: f64) -> Result<(), LambertError> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(LambertError::NonFiniteInput { parameter, value })
    }
}

fn validate_finite_vector(prefix: &'static str, value: Vec3) -> Result<(), LambertError> {
    validate_finite_scalar(component_name(prefix, "x"), value[0])?;
    validate_finite_scalar(component_name(prefix, "y"), value[1])?;
    validate_finite_scalar(component_name(prefix, "z"), value[2])
}

fn component_name(prefix: &'static str, component: &'static str) -> &'static str {
    match (prefix, component) {
        ("r1_km", "x") => "r1_km.x",
        ("r1_km", "y") => "r1_km.y",
        ("r1_km", "z") => "r1_km.z",
        ("r2_km", "x") => "r2_km.x",
        ("r2_km", "y") => "r2_km.y",
        ("r2_km", "z") => "r2_km.z",
        _ => prefix,
    }
}
