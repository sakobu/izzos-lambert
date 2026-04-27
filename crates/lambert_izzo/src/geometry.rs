//! Lambert problem geometry: chord, semi-perimeter, λ, and transfer-plane basis.
//!
//! Constructed once per [`crate::lambert`] call and threaded through the
//! root-finding and velocity-reconstruction stages — eliminates re-derivation
//! drift across modules.

// See `vec3.rs` for the rationale on this `allow(unused_imports)`.
#[allow(unused_imports)]
use num_traits::Float as _;

use crate::constants::{COLINEARITY_TOL, MIN_POSITION_NORM};
use crate::error::{LambertError, NonFiniteParameter, Position};
use crate::vec3::{self, Vec3};

/// Pre-computed geometry for a Lambert boundary problem.
///
/// All scalars and unit vectors derive from `(r1, r2, tof, mu, way)`. The
/// solver kernels consume `lambda` and `big_t` (non-dimensional TOF); the
/// velocity reconstruction in [`crate::lambert`] consumes the rest.
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
    /// - [`LambertError::NonPositiveTimeOfFlight`] — `tof <= 0`.
    /// - [`LambertError::NonPositiveMu`] — `mu <= 0`.
    /// - [`LambertError::DegeneratePositionVector`] — `|r1|` or `|r2|`
    ///   below [`MIN_POSITION_NORM`].
    /// - [`LambertError::CollinearGeometry`] — `|r1 × r2| / (|r1| · |r2|)`
    ///   below [`COLINEARITY_TOL`].
    #[allow(clippy::similar_names)] // ir1/ir2 are radial unit vectors, it1/it2 tangential — Izzo Eq. 5–7.
    pub(crate) fn from_inputs(
        r1: Vec3,
        r2: Vec3,
        tof: f64,
        mu: f64,
        way: crate::TransferWay,
    ) -> Result<Self, LambertError> {
        validate_finite_r1(r1)?;
        validate_finite_r2(r2)?;
        validate_finite_scalar(NonFiniteParameter::Tof, tof)?;
        validate_finite_scalar(NonFiniteParameter::Mu, mu)?;

        if tof <= 0.0 {
            return Err(LambertError::NonPositiveTimeOfFlight { tof });
        }
        if mu <= 0.0 {
            return Err(LambertError::NonPositiveMu { mu });
        }

        let r1n = vec3::norm(r1);
        let r2n = vec3::norm(r2);
        if r1n < MIN_POSITION_NORM {
            return Err(LambertError::DegeneratePositionVector {
                position: Position::R1,
                norm: r1n,
            });
        }
        if r2n < MIN_POSITION_NORM {
            return Err(LambertError::DegeneratePositionVector {
                position: Position::R2,
                norm: r2n,
            });
        }

        let chord = vec3::sub(r2, r1);
        let c = vec3::norm(chord);
        let s = 0.5 * (r1n + r2n + c);

        let ir1 = vec3::scale(r1, 1.0 / r1n);
        let ir2 = vec3::scale(r2, 1.0 / r2n);
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

        let big_t = (2.0 * mu / (s * s * s)).sqrt() * tof;
        let gamma = (mu * s / 2.0).sqrt();
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

fn validate_finite_scalar(parameter: NonFiniteParameter, value: f64) -> Result<(), LambertError> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(LambertError::NonFiniteInput { parameter, value })
    }
}

fn validate_finite_r1(value: Vec3) -> Result<(), LambertError> {
    validate_finite_scalar(NonFiniteParameter::R1X, value[0])?;
    validate_finite_scalar(NonFiniteParameter::R1Y, value[1])?;
    validate_finite_scalar(NonFiniteParameter::R1Z, value[2])
}

fn validate_finite_r2(value: Vec3) -> Result<(), LambertError> {
    validate_finite_scalar(NonFiniteParameter::R2X, value[0])?;
    validate_finite_scalar(NonFiniteParameter::R2Y, value[1])?;
    validate_finite_scalar(NonFiniteParameter::R2Z, value[2])
}
