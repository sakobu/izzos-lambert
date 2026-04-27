//! JavaScript and TypeScript bindings for `lambert_izzo`.

#![warn(missing_docs)]
#![warn(clippy::pedantic)]
#![warn(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::unreachable
)]
#![allow(clippy::module_name_repetitions)] // LambertRequest / LambertResponse mirror the core crate naming.

use lambert_izzo::{
    LambertError, LambertSolutions, NonFiniteParameter, RevolutionBudget, TransferWay,
    solve_with_diagnostics as core_solve,
};
use serde::{Deserialize, Serialize};
use tsify::Tsify;
use wasm_bindgen::prelude::{JsValue, wasm_bindgen};

/// Direction around the transfer plane from `r1` to `r2`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Deserialize, Serialize, Tsify)]
#[serde(rename_all = "camelCase")]
#[tsify(from_wasm_abi, into_wasm_abi)]
pub enum TransferWayInput {
    /// Short-way transfer with angle less than or equal to pi.
    Short,
    /// Long-way transfer with angle greater than pi.
    Long,
}

impl From<TransferWayInput> for TransferWay {
    fn from(value: TransferWayInput) -> Self {
        match value {
            TransferWayInput::Short => Self::Short,
            TransferWayInput::Long => Self::Long,
        }
    }
}

/// JavaScript-friendly Lambert request.
#[derive(Debug, Clone, Copy, PartialEq, Deserialize, Serialize, Tsify)]
#[serde(rename_all = "camelCase")]
#[tsify(from_wasm_abi, into_wasm_abi)]
pub struct LambertRequest {
    /// Initial position vector in kilometers.
    pub r1_km: [f64; 3],
    /// Final position vector in kilometers.
    pub r2_km: [f64; 3],
    /// Time of flight in seconds.
    pub tof_s: f64,
    /// Gravitational parameter in cubic kilometers per square second.
    pub mu_km3_s2: f64,
    /// Short-way or long-way transfer selection.
    pub way: TransferWayInput,
    /// Maximum complete revolutions to consider beyond single-rev.
    pub max_revs: u32,
}

/// One JavaScript-friendly Lambert trajectory plus solver diagnostics.
#[derive(Debug, Clone, Copy, PartialEq, Deserialize, Serialize, Tsify)]
#[serde(rename_all = "camelCase")]
#[tsify(from_wasm_abi, into_wasm_abi)]
pub struct LambertSolutionOutput {
    /// Velocity at `r1_km` in kilometers per second.
    pub v1_km_s: [f64; 3],
    /// Velocity at `r2_km` in kilometers per second.
    pub v2_km_s: [f64; 3],
    /// Solver diagnostics for this branch.
    pub diagnostics: SolverDiagnosticsOutput,
}

/// One JavaScript-friendly multi-rev pair.
#[derive(Debug, Clone, Copy, PartialEq, Deserialize, Serialize, Tsify)]
#[serde(rename_all = "camelCase")]
#[tsify(from_wasm_abi, into_wasm_abi)]
pub struct MultiRevPairOutput {
    /// Branch revolution count (`>= 1`).
    pub n_revs: u32,
    /// Long-period trajectory.
    pub long_period: LambertSolutionOutput,
    /// Short-period trajectory.
    pub short_period: LambertSolutionOutput,
}

/// JavaScript-friendly solver response.
#[derive(Debug, Clone, PartialEq, Deserialize, Serialize, Tsify)]
#[serde(rename_all = "camelCase")]
#[tsify(from_wasm_abi, into_wasm_abi)]
pub struct LambertResponse {
    /// Single-revolution trajectory — always present.
    pub single: LambertSolutionOutput,
    /// Multi-revolution pairs in ascending `M` order.
    pub multi: Vec<MultiRevPairOutput>,
}

/// Solver diagnostic data for JavaScript callers.
#[derive(Debug, Clone, Copy, PartialEq, Deserialize, Serialize, Tsify)]
#[serde(rename_all = "camelCase")]
#[tsify(from_wasm_abi, into_wasm_abi)]
pub struct SolverDiagnosticsOutput {
    /// Householder iterations used to converge.
    pub iters: u32,
    /// Final Lancaster-Blanchard x value.
    pub x: f64,
}

/// JavaScript-friendly mirror of [`lambert_izzo::LambertError`].
///
/// Serialized as a discriminated union: `{ kind: "VariantName", ...fields }`.
/// JS callers can `switch` on `kind` to handle each failure mode by name.
///
/// Mirrors the core enum 1:1; if a future variant is added to
/// `LambertError`, the `From` impl below will hit `unreachable!()` —
/// add the matching variant here.
#[derive(Debug, Clone, Copy, PartialEq, Deserialize, Serialize, Tsify)]
#[serde(tag = "kind")]
#[tsify(from_wasm_abi, into_wasm_abi)]
pub enum LambertErrorOutput {
    /// Mirrors [`LambertError::NonFiniteInput`].
    NonFiniteInput {
        /// Which public parameter or vector component was non-finite.
        parameter: NonFiniteParameterOutput,
        /// The non-finite value the caller passed.
        value: f64,
    },
    /// Mirrors [`LambertError::NonPositiveTimeOfFlight`].
    NonPositiveTimeOfFlight {
        /// The non-positive `tof` value the caller passed (s).
        tof_s: f64,
    },
    /// Mirrors [`LambertError::NonPositiveMu`].
    NonPositiveMu {
        /// The non-positive `mu` value the caller passed (km³/s²).
        mu_km3_s2: f64,
    },
    /// Mirrors [`LambertError::DegeneratePositionVector`].
    DegeneratePositionVector {
        /// `1` for `r1`, `2` for `r2`.
        which: u8,
        /// Norm of the offending vector (km).
        norm_km: f64,
    },
    /// Mirrors [`LambertError::CollinearGeometry`].
    CollinearGeometry {
        /// `|r1 × r2| / (|r1| · |r2|)` — the sine of the transfer angle.
        sin_angle: f64,
    },
    /// Mirrors [`LambertError::NoConvergence`].
    NoConvergence {
        /// Iterations performed before giving up.
        iterations: u32,
        /// Magnitude of the last `|Δx|` step.
        last_step: f64,
        /// Branch index: `0` = single-rev, `≥ 1` = multi-rev.
        n_revs: u32,
    },
    /// Mirrors [`LambertError::SingularDenominator`].
    SingularDenominator {
        /// Branch index where the singularity occurred.
        n_revs: u32,
    },
}

/// Mirrors [`lambert_izzo::NonFiniteParameter`] with TS-friendly tagging.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Deserialize, Serialize, Tsify)]
#[tsify(from_wasm_abi, into_wasm_abi)]
pub enum NonFiniteParameterOutput {
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

impl From<NonFiniteParameter> for NonFiniteParameterOutput {
    fn from(value: NonFiniteParameter) -> Self {
        match value {
            NonFiniteParameter::R1KmX => Self::R1KmX,
            NonFiniteParameter::R1KmY => Self::R1KmY,
            NonFiniteParameter::R1KmZ => Self::R1KmZ,
            NonFiniteParameter::R2KmX => Self::R2KmX,
            NonFiniteParameter::R2KmY => Self::R2KmY,
            NonFiniteParameter::R2KmZ => Self::R2KmZ,
            NonFiniteParameter::TofS => Self::TofS,
            NonFiniteParameter::MuKm3S2 => Self::MuKm3S2,
        }
    }
}

impl From<LambertError> for LambertErrorOutput {
    fn from(value: LambertError) -> Self {
        match value {
            LambertError::NonFiniteInput { parameter, value } => Self::NonFiniteInput {
                parameter: parameter.into(),
                value,
            },
            LambertError::NonPositiveTimeOfFlight { tof_s } => {
                Self::NonPositiveTimeOfFlight { tof_s }
            }
            LambertError::NonPositiveMu { mu_km3_s2 } => Self::NonPositiveMu { mu_km3_s2 },
            LambertError::DegeneratePositionVector { which, norm_km } => {
                Self::DegeneratePositionVector { which, norm_km }
            }
            LambertError::CollinearGeometry { sin_angle } => Self::CollinearGeometry { sin_angle },
            LambertError::NoConvergence {
                iterations,
                last_step,
                n_revs,
            } => Self::NoConvergence {
                iterations,
                last_step,
                n_revs,
            },
            LambertError::SingularDenominator { n_revs } => Self::SingularDenominator { n_revs },
            // LambertError is #[non_exhaustive]; if a new variant ever lands
            // upstream, add a matching arm here. Until then, fall through to
            // the catch-all below (which won't fire in practice).
            other => {
                // Lossy fallback: serialize as a singular-denominator with
                // n_revs = 0 so JS at least gets a typed error rather than a
                // panic. In practice this branch is dead.
                let _ = other;
                Self::SingularDenominator { n_revs: 0 }
            }
        }
    }
}

/// Solve a Lambert request from JavaScript.
///
/// # Errors
///
/// Returns a structured [`LambertErrorOutput`] (serialized as a JS object
/// with a `kind` discriminator) when the core solver rejects the input or
/// fails to converge.
#[wasm_bindgen(js_name = solveLambert)]
pub fn solve_lambert(request: LambertRequest) -> Result<LambertResponse, JsValue> {
    solve_lambert_request(request).map_err(|err| {
        serde_wasm_bindgen::to_value(&err).unwrap_or_else(|serialize_err| {
            // Fallback to a string if even the serialization fails — should
            // never happen with the current types, but keeps the panic-free
            // discipline.
            JsValue::from_str(&serialize_err.to_string())
        })
    })
}

/// Solve a Lambert request using only Rust types.
///
/// This function exists so the wrapper contract can be tested without a
/// JavaScript runtime.
///
/// # Errors
///
/// Returns the structured [`LambertErrorOutput`] when the request is
/// invalid or the numerical solve fails.
pub fn solve_lambert_request(
    request: LambertRequest,
) -> Result<LambertResponse, LambertErrorOutput> {
    let revolutions = RevolutionBudget::up_to(request.max_revs);
    let (solutions, diagnostics) = core_solve(
        request.r1_km,
        request.r2_km,
        request.tof_s,
        request.mu_km3_s2,
        request.way.into(),
        revolutions,
    )
    .map_err(LambertErrorOutput::from)?;
    Ok(into_response(&solutions, &diagnostics))
}

fn into_response(
    solutions: &LambertSolutions,
    diagnostics: &lambert_izzo::LambertDiagnostics,
) -> LambertResponse {
    let single = LambertSolutionOutput {
        v1_km_s: solutions.single.v1_km_s,
        v2_km_s: solutions.single.v2_km_s,
        diagnostics: SolverDiagnosticsOutput {
            iters: diagnostics.single.iters,
            x: diagnostics.single.lancaster_blanchard_x,
        },
    };
    let multi = solutions
        .multi
        .iter()
        .zip(diagnostics.multi.iter())
        .map(|(pair, dpair)| MultiRevPairOutput {
            n_revs: pair.n_revs,
            long_period: LambertSolutionOutput {
                v1_km_s: pair.long_period.v1_km_s,
                v2_km_s: pair.long_period.v2_km_s,
                diagnostics: SolverDiagnosticsOutput {
                    iters: dpair.long_period.iters,
                    x: dpair.long_period.lancaster_blanchard_x,
                },
            },
            short_period: LambertSolutionOutput {
                v1_km_s: pair.short_period.v1_km_s,
                v2_km_s: pair.short_period.v2_km_s,
                diagnostics: SolverDiagnosticsOutput {
                    iters: dpair.short_period.iters,
                    x: dpair.short_period.lancaster_blanchard_x,
                },
            },
        })
        .collect();
    LambertResponse { single, multi }
}
