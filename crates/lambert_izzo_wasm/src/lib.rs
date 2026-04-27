//! JavaScript and TypeScript bindings for `lambert_izzo`.

#![warn(missing_docs)]
#![warn(clippy::pedantic)]
#![warn(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::unreachable
)]
#![allow(clippy::module_name_repetitions)] // LambertRequest and LambertResponse mirror the core crate naming.

use lambert_izzo::{RevolutionBudget, TransferWay, lambert};
use nalgebra::Vector3;
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

/// JavaScript-friendly solver response.
#[derive(Debug, Clone, PartialEq, Deserialize, Serialize, Tsify)]
#[serde(rename_all = "camelCase")]
#[tsify(from_wasm_abi, into_wasm_abi)]
pub struct LambertResponse {
    /// All reachable Lambert solutions.
    pub solutions: Vec<LambertSolutionOutput>,
}

/// One JavaScript-friendly Lambert solution.
#[derive(Debug, Clone, PartialEq, Deserialize, Serialize, Tsify)]
#[serde(rename_all = "camelCase")]
#[tsify(from_wasm_abi, into_wasm_abi)]
pub struct LambertSolutionOutput {
    /// Velocity at `r1_km` in kilometers per second.
    pub v1_km_s: [f64; 3],
    /// Velocity at `r2_km` in kilometers per second.
    pub v2_km_s: [f64; 3],
    /// Number of complete revolutions for this branch.
    pub n_revs: u32,
    /// Solver diagnostics for this branch.
    pub diagnostics: SolverDiagnosticsOutput,
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

/// Solve a Lambert request from JavaScript.
///
/// # Errors
///
/// Returns a JavaScript string when the core solver rejects the input or fails
/// to converge.
#[wasm_bindgen(js_name = solveLambert)]
pub fn solve_lambert(request: LambertRequest) -> Result<LambertResponse, JsValue> {
    solve_lambert_request(request).map_err(|error| JsValue::from_str(&error))
}

/// Solve a Lambert request using only Rust types.
///
/// This function exists so the wrapper contract can be tested without a
/// JavaScript runtime.
///
/// # Errors
///
/// Returns the core solver error message when the request is invalid or the
/// numerical solve fails.
pub fn solve_lambert_request(request: LambertRequest) -> Result<LambertResponse, String> {
    let r1_km = array_to_vector3(request.r1_km);
    let r2_km = array_to_vector3(request.r2_km);
    let revolutions = RevolutionBudget::up_to(request.max_revs);
    let solutions = lambert(
        r1_km,
        r2_km,
        request.tof_s,
        request.mu_km3_s2,
        request.way.into(),
        revolutions,
    )
    .map_err(|error| error.to_string())?
    .into_iter()
    .map(|solution| LambertSolutionOutput {
        v1_km_s: vector3_to_array(solution.v1_km_s),
        v2_km_s: vector3_to_array(solution.v2_km_s),
        n_revs: solution.n_revs,
        diagnostics: SolverDiagnosticsOutput {
            iters: solution.diagnostics.iters,
            x: solution.diagnostics.lancaster_blanchard_x,
        },
    })
    .collect();

    Ok(LambertResponse { solutions })
}

fn array_to_vector3(value: [f64; 3]) -> Vector3<f64> {
    Vector3::new(value[0], value[1], value[2])
}

fn vector3_to_array(value: Vector3<f64>) -> [f64; 3] {
    [value.x, value.y, value.z]
}
