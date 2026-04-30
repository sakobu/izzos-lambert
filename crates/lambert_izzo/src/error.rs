//! Error type returned by [`crate::lambert`].

/// Identifies which of the two endpoint positions an error refers to.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum Position {
    /// The initial position, `r1`.
    R1,
    /// The final position, `r2`.
    R2,
}

impl core::fmt::Display for Position {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.write_str(match self {
            Self::R1 => "r1",
            Self::R2 => "r2",
        })
    }
}

/// Identifies which public input was rejected by the finiteness check.
///
/// A typed alternative to a free-form `&'static str`, so the error type
/// stays `Deserialize`-friendly under the `serde` feature without falling
/// back to allocated strings.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum NonFiniteParameter {
    /// `r1.x`
    R1X,
    /// `r1.y`
    R1Y,
    /// `r1.z`
    R1Z,
    /// `r2.x`
    R2X,
    /// `r2.y`
    R2Y,
    /// `r2.z`
    R2Z,
    /// `tof`
    Tof,
    /// `mu`
    Mu,
}

impl NonFiniteParameter {
    /// Public-API name of the offending parameter.
    #[must_use]
    pub fn as_str(self) -> &'static str {
        match self {
            Self::R1X => "r1.x",
            Self::R1Y => "r1.y",
            Self::R1Z => "r1.z",
            Self::R2X => "r2.x",
            Self::R2Y => "r2.y",
            Self::R2Z => "r2.z",
            Self::Tof => "tof",
            Self::Mu => "mu",
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
/// All scalar fields use the crate's documented unit convention:
/// positions in km, time of flight in seconds, gravitational parameter
/// in km³/s². Unitless fields (`sin_angle`, `last_step` —
/// Izzo's dimensionless `x`-step) carry no unit at all.
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
    #[error("non-positive time of flight: tof = {tof}")]
    NonPositiveTimeOfFlight {
        /// The non-positive `tof` value the caller passed (s).
        tof: f64,
    },

    /// Gravitational parameter must be strictly positive.
    #[error("non-positive gravitational parameter: mu = {mu}")]
    NonPositiveMu {
        /// The non-positive `mu` value the caller passed (km³/s²).
        mu: f64,
    },

    /// One position vector has near-zero norm; geometry undefined.
    ///
    /// Triggered when `|r|` is below an internal floor (`1e-15` km).
    #[error("degenerate position vector {position}: norm = {norm}")]
    DegeneratePositionVector {
        /// Which of the two positions is degenerate.
        position: Position,
        /// Norm of the offending vector (km).
        norm: f64,
    },

    /// `r1` and `r2` are colinear; the transfer plane is undefined.
    ///
    /// Triggered when `|r1 × r2| / (|r1| · |r2|)` is below an internal
    /// floor (`1e-15`, scale-invariant).
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
