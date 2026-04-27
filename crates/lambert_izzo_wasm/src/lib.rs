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
    LambertError, LambertSolutions, NonFiniteParameter, Position, RevolutionBudget, TransferWay,
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

/// JavaScript-friendly Lambert request. All scalars use the crate's
/// documented km / s / km³/s² convention.
#[derive(Debug, Clone, Copy, PartialEq, Deserialize, Serialize, Tsify)]
#[serde(rename_all = "camelCase")]
#[tsify(from_wasm_abi, into_wasm_abi)]
pub struct LambertRequest {
    /// Initial position vector.
    pub r1: [f64; 3],
    /// Final position vector.
    pub r2: [f64; 3],
    /// Time of flight.
    pub tof: f64,
    /// Gravitational parameter of the central body.
    pub mu: f64,
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
    /// Velocity at `r1`.
    pub v1: [f64; 3],
    /// Velocity at `r2`.
    pub v2: [f64; 3],
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

/// Mirrors [`lambert_izzo::Position`] with TS-friendly tagging.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Deserialize, Serialize, Tsify)]
#[tsify(from_wasm_abi, into_wasm_abi)]
pub enum PositionOutput {
    /// The initial position, `r1`.
    R1,
    /// The final position, `r2`.
    R2,
}

impl From<Position> for PositionOutput {
    fn from(value: Position) -> Self {
        match value {
            Position::R1 => Self::R1,
            Position::R2 => Self::R2,
        }
    }
}

/// Mirrors [`lambert_izzo::NonFiniteParameter`] with TS-friendly tagging.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Deserialize, Serialize, Tsify)]
#[tsify(from_wasm_abi, into_wasm_abi)]
pub enum NonFiniteParameterOutput {
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

impl From<NonFiniteParameter> for NonFiniteParameterOutput {
    fn from(value: NonFiniteParameter) -> Self {
        match value {
            NonFiniteParameter::R1X => Self::R1X,
            NonFiniteParameter::R1Y => Self::R1Y,
            NonFiniteParameter::R1Z => Self::R1Z,
            NonFiniteParameter::R2X => Self::R2X,
            NonFiniteParameter::R2Y => Self::R2Y,
            NonFiniteParameter::R2Z => Self::R2Z,
            NonFiniteParameter::Tof => Self::Tof,
            NonFiniteParameter::Mu => Self::Mu,
        }
    }
}

/// JavaScript-friendly mirror of [`lambert_izzo::LambertError`].
///
/// Serialized as a discriminated union: `{ kind: "VariantName", ...fields }`.
/// JS callers can `switch` on `kind` to handle each failure mode by name.
///
/// Mirrors the core enum 1:1; if a future variant is added to
/// `LambertError`, the `From` impl below will fall through to the catch-
/// all (`SingularDenominator { n_revs: 0 }`) — add the matching variant
/// here when that happens.
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
        /// The non-positive `tof` value.
        tof: f64,
    },
    /// Mirrors [`LambertError::NonPositiveMu`].
    NonPositiveMu {
        /// The non-positive `mu` value.
        mu: f64,
    },
    /// Mirrors [`LambertError::DegeneratePositionVector`].
    DegeneratePositionVector {
        /// Which of the two positions is degenerate.
        position: PositionOutput,
        /// Norm of the offending vector.
        norm: f64,
    },
    /// Mirrors [`LambertError::CollinearGeometry`].
    CollinearGeometry {
        /// `|r1 × r2| / (|r1| · |r2|)` — sine of the transfer angle.
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

impl From<LambertError> for LambertErrorOutput {
    fn from(value: LambertError) -> Self {
        match value {
            LambertError::NonFiniteInput { parameter, value } => Self::NonFiniteInput {
                parameter: parameter.into(),
                value,
            },
            LambertError::NonPositiveTimeOfFlight { tof } => Self::NonPositiveTimeOfFlight { tof },
            LambertError::NonPositiveMu { mu } => Self::NonPositiveMu { mu },
            LambertError::DegeneratePositionVector { position, norm } => {
                Self::DegeneratePositionVector {
                    position: position.into(),
                    norm,
                }
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
            // a catch-all that won't fire in practice.
            other => {
                let _ = other;
                Self::SingularDenominator { n_revs: 0 }
            }
        }
    }
}

/// Solve a Lambert request from JavaScript.
///
/// # Invariants
///
/// Same as [`lambert_izzo::lambert`]. Inputs come through
/// [`LambertRequest`]; `tof > 0`, `mu > 0`, finite vectors, transfer
/// angle ∉ {0, π}.
///
/// # Validity / near-degenerate behavior
///
/// Same as [`lambert_izzo::lambert`].
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
/// # Invariants
///
/// Same as [`lambert_izzo::lambert`]. Inputs come through
/// [`LambertRequest`]; `tof > 0`, `mu > 0`, finite vectors, transfer
/// angle ∉ {0, π}.
///
/// # Validity / near-degenerate behavior
///
/// Same as [`lambert_izzo::lambert`].
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
        request.r1,
        request.r2,
        request.tof,
        request.mu,
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
        v1: solutions.single.v1,
        v2: solutions.single.v2,
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
                v1: pair.long_period.v1,
                v2: pair.long_period.v2,
                diagnostics: SolverDiagnosticsOutput {
                    iters: dpair.long_period.iters,
                    x: dpair.long_period.lancaster_blanchard_x,
                },
            },
            short_period: LambertSolutionOutput {
                v1: pair.short_period.v1,
                v2: pair.short_period.v2,
                diagnostics: SolverDiagnosticsOutput {
                    iters: dpair.short_period.iters,
                    x: dpair.short_period.lancaster_blanchard_x,
                },
            },
        })
        .collect();
    LambertResponse { single, multi }
}
