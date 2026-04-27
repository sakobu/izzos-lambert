//! Error type returned by [`crate::lambert`].

/// Identifies which public input was rejected by the finiteness check.
///
/// A typed alternative to a free-form `&'static str`, so the error type
/// stays `Deserialize`-friendly under the `serde` feature without falling
/// back to allocated strings.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum NonFiniteParameter {
    /// `r1_km.x`
    R1KmX,
    /// `r1_km.y`
    R1KmY,
    /// `r1_km.z`
    R1KmZ,
    /// `r2_km.x`
    R2KmX,
    /// `r2_km.y`
    R2KmY,
    /// `r2_km.z`
    R2KmZ,
    /// `tof_s`
    TofS,
    /// `mu_km3_s2`
    MuKm3S2,
}

impl NonFiniteParameter {
    /// String name of the offending parameter, matching the public API name.
    #[must_use]
    pub fn as_str(self) -> &'static str {
        match self {
            Self::R1KmX => "r1_km.x",
            Self::R1KmY => "r1_km.y",
            Self::R1KmZ => "r1_km.z",
            Self::R2KmX => "r2_km.x",
            Self::R2KmY => "r2_km.y",
            Self::R2KmZ => "r2_km.z",
            Self::TofS => "tof_s",
            Self::MuKm3S2 => "mu_km3_s2",
        }
    }
}

impl core::fmt::Display for NonFiniteParameter {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.write_str(self.as_str())
    }
}

/// Failure modes of the Izzo Lambert solver.
///
/// Field units follow the crate's SI convention: `_km` for lengths, `_s` for
/// times, `_km3_s2` for the gravitational parameter. Unitless fields
/// (`sin_angle`, `last_step` — Izzo's dimensionless `x`-step) carry no suffix.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, thiserror::Error)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum LambertError {
    /// One public input was `NaN`, `+inf`, or `-inf`.
    #[error("non-finite input: {parameter} = {value}")]
    NonFiniteInput {
        /// Which public parameter or vector component was non-finite.
        parameter: NonFiniteParameter,
        /// The non-finite value the caller passed.
        value: f64,
    },

    /// Time of flight must be strictly positive.
    #[error("non-positive time of flight: tof_s = {tof_s}")]
    NonPositiveTimeOfFlight {
        /// The non-positive `tof` value the caller passed (s).
        tof_s: f64,
    },

    /// Gravitational parameter must be strictly positive.
    #[error("non-positive gravitational parameter: mu_km3_s2 = {mu_km3_s2}")]
    NonPositiveMu {
        /// The non-positive `mu` value the caller passed (km³/s²).
        mu_km3_s2: f64,
    },

    /// One position vector has near-zero norm; geometry undefined.
    ///
    /// Triggered when `|r_which|` is below
    /// [`crate::constants::MIN_POSITION_NORM_KM`].
    #[error("degenerate position vector r{which}: norm_km = {norm_km}")]
    DegeneratePositionVector {
        /// `1` for `r1`, `2` for `r2`.
        which: u8,
        /// Norm of the offending vector (km).
        norm_km: f64,
    },

    /// `r1` and `r2` are colinear; the transfer plane is undefined.
    ///
    /// Triggered when `|r1 × r2| / (|r1| · |r2|)` is below
    /// [`crate::constants::COLINEARITY_TOL`].
    #[error("colinear position vectors: |r1 × r2| / (|r1| |r2|) = {sin_angle:.3e}")]
    CollinearGeometry {
        /// `|r1 × r2| / (|r1| · |r2|)` — the sine of the transfer angle (unitless).
        ///
        /// Stored rather than the angle itself because `asin` near `0` or `π`
        /// is the noisier of the two — the sine is what the check used.
        sin_angle: f64,
    },

    /// Householder iteration did not reach the configured tolerance.
    #[error(
        "Householder did not converge after {iterations} iters \
         (last |Δx| = {last_step:.3e}, branch M = {n_revs})"
    )]
    NoConvergence {
        /// Iterations performed before giving up.
        iterations: u32,
        /// Magnitude of the last `|Δx|` step (Izzo's `x`-space, unitless).
        last_step: f64,
        /// Branch index: `0` = single-rev, `≥ 1` = multi-rev.
        n_revs: u32,
    },

    /// Householder denominator collapsed to zero — algebraic singularity,
    /// distinct from slow iterative convergence.
    #[error("Householder denominator vanished on branch M = {n_revs}")]
    SingularDenominator {
        /// Branch index where the singularity occurred.
        n_revs: u32,
    },
}
