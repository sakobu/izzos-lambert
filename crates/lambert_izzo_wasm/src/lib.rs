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
    lambert as core_solve,
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

/// One JavaScript-friendly Lambert trajectory.
///
/// Diagnostics are no longer attached per-solution — see
/// [`LambertResponse::diagnostics`] for the per-branch iteration counts.
#[derive(Debug, Clone, Copy, PartialEq, Deserialize, Serialize, Tsify)]
#[serde(rename_all = "camelCase")]
#[tsify(from_wasm_abi, into_wasm_abi)]
pub struct LambertSolutionOutput {
    /// Velocity at `r1`.
    pub v1: [f64; 3],
    /// Velocity at `r2`.
    pub v2: [f64; 3],
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
///
/// Always carries the single-revolution trajectory, every reachable
/// multi-rev pair (in ascending `M` order), and the per-branch
/// [`LambertDiagnosticsOutput`] mirroring the structure 1:1.
#[derive(Debug, Clone, PartialEq, Deserialize, Serialize, Tsify)]
#[serde(rename_all = "camelCase")]
#[tsify(from_wasm_abi, into_wasm_abi)]
pub struct LambertResponse {
    /// Single-revolution trajectory — always present.
    pub single: LambertSolutionOutput,
    /// Multi-revolution pairs in ascending `M` order.
    pub multi: Vec<MultiRevPairOutput>,
    /// Per-branch solver diagnostics, structurally parallel to
    /// `single` / `multi`.
    pub diagnostics: LambertDiagnosticsOutput,
}

/// Solver diagnostic data for one converged branch.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Deserialize, Serialize, Tsify)]
#[serde(rename_all = "camelCase")]
#[tsify(from_wasm_abi, into_wasm_abi)]
pub struct SolverDiagnosticsOutput {
    /// Householder iterations used to converge.
    pub iters: u32,
}

/// Diagnostics for one multi-rev pair.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Deserialize, Serialize, Tsify)]
#[serde(rename_all = "camelCase")]
#[tsify(from_wasm_abi, into_wasm_abi)]
pub struct MultiRevPairDiagnosticsOutput {
    /// Branch revolution count (`>= 1`).
    pub n_revs: u32,
    /// Long-period branch diagnostics.
    pub long_period: SolverDiagnosticsOutput,
    /// Short-period branch diagnostics.
    pub short_period: SolverDiagnosticsOutput,
}

/// Diagnostics structure mirroring [`LambertResponse`].
#[derive(Debug, Clone, PartialEq, Deserialize, Serialize, Tsify)]
#[serde(rename_all = "camelCase")]
#[tsify(from_wasm_abi, into_wasm_abi)]
pub struct LambertDiagnosticsOutput {
    /// Single-rev solver diagnostics.
    pub single: SolverDiagnosticsOutput,
    /// Multi-rev pair diagnostics in the same order as
    /// [`LambertResponse::multi`].
    pub multi: Vec<MultiRevPairDiagnosticsOutput>,
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
/// `LambertError` is `#[non_exhaustive]`, so a literal exhaustive match is
/// not available to downstream crates. The [`Unknown`] variant carries the
/// upstream error's `Display` rendering for any future variant we have not
/// yet mirrored — JS callers should treat it as fatal and surface the
/// message directly rather than silently mapping it to a wrong `kind`.
///
/// [`Unknown`]: LambertErrorOutput::Unknown
#[derive(Debug, Clone, PartialEq, Deserialize, Serialize, Tsify)]
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
    /// A new [`LambertError`] variant has been added upstream that this
    /// wasm adapter does not yet mirror. The `message` field is the
    /// upstream error's `Display` rendering; JS callers should report it
    /// verbatim. When this fires, add a matching arm to
    /// `From<LambertError> for LambertErrorOutput`.
    Unknown {
        /// `Display` text of the upstream `LambertError`.
        message: String,
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
            // LambertError is #[non_exhaustive]; new variants land in
            // `Unknown` with the upstream Display text rather than being
            // silently mislabeled. Add a matching arm above when extending.
            other => Self::Unknown {
                message: format!("{other}"),
            },
        }
    }
}

/// One element of a batch response: success carries the response, failure
/// carries the structured error.
///
/// Serialized as a tagged union: `{ kind: "ok", response: ... }` or
/// `{ kind: "err", error: ... }`. JS can narrow on `result.kind`.
#[derive(Debug, Clone, PartialEq, Deserialize, Serialize, Tsify)]
#[serde(tag = "kind", rename_all = "lowercase")]
#[tsify(from_wasm_abi, into_wasm_abi)]
pub enum BatchResult {
    /// The solver succeeded; `response` carries the trajectories and
    /// diagnostics.
    Ok {
        /// Successful response.
        response: LambertResponse,
    },
    /// The solver rejected the input or failed to converge; `error`
    /// carries the structured failure.
    Err {
        /// Structured failure.
        error: LambertErrorOutput,
    },
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

/// Solve a batch of Lambert requests from JavaScript.
///
/// Returns one [`BatchResult`] per request, in input order. A failed
/// request does not stop processing of the rest — each element is
/// independent.
///
/// # Invariants / Validity / Errors
///
/// Per-element, same as [`solve_lambert`].
#[wasm_bindgen(js_name = solveLambertBatch)]
#[must_use]
pub fn solve_lambert_batch(requests: Vec<LambertRequest>) -> Vec<BatchResult> {
    requests
        .into_iter()
        .map(|request| match solve_lambert_request(request) {
            Ok(response) => BatchResult::Ok { response },
            Err(error) => BatchResult::Err { error },
        })
        .collect()
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
    let solutions = core_solve(&lambert_izzo::LambertInput {
        r1: request.r1,
        r2: request.r2,
        tof: request.tof,
        mu: request.mu,
        way: request.way.into(),
        revolutions,
    })
    .map_err(LambertErrorOutput::from)?;
    Ok(into_response(&solutions))
}

fn into_response(solutions: &LambertSolutions) -> LambertResponse {
    let single = LambertSolutionOutput {
        v1: solutions.single.v1,
        v2: solutions.single.v2,
    };

    let n = solutions.multi.len();
    let mut multi = Vec::with_capacity(n);
    let mut multi_diagnostics = Vec::with_capacity(n);
    for (pair, dpair) in solutions
        .multi
        .iter()
        .zip(solutions.diagnostics.multi.iter())
    {
        multi.push(MultiRevPairOutput {
            n_revs: pair.n_revs,
            long_period: LambertSolutionOutput {
                v1: pair.long_period.v1,
                v2: pair.long_period.v2,
            },
            short_period: LambertSolutionOutput {
                v1: pair.short_period.v1,
                v2: pair.short_period.v2,
            },
        });
        multi_diagnostics.push(MultiRevPairDiagnosticsOutput {
            n_revs: dpair.n_revs,
            long_period: SolverDiagnosticsOutput {
                iters: dpair.long_period.iters,
            },
            short_period: SolverDiagnosticsOutput {
                iters: dpair.short_period.iters,
            },
        });
    }

    LambertResponse {
        single,
        multi,
        diagnostics: LambertDiagnosticsOutput {
            single: SolverDiagnosticsOutput {
                iters: solutions.diagnostics.single.iters,
            },
            multi: multi_diagnostics,
        },
    }
}
