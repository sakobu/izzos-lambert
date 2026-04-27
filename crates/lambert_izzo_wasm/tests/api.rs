use lambert_izzo_test_support::bodies::MU_EARTH;
use lambert_izzo_wasm::{
    LambertErrorOutput, LambertRequest, NonFiniteParameterOutput, PositionOutput,
    TransferWayInput, solve_lambert_request,
};

#[test]
fn single_rev_request_returns_js_friendly_vectors() {
    let request = LambertRequest {
        r1: [7000.0, 0.0, 0.0],
        r2: [0.0, 7000.0, 0.0],
        tof: core::f64::consts::PI / 2.0 * (7000.0_f64.powi(3) / MU_EARTH).sqrt(),
        mu: MU_EARTH,
        way: TransferWayInput::Short,
        max_revs: 0,
    };

    let response = solve_lambert_request(request).unwrap();

    assert!(response.multi.is_empty());
    let single = response.single;
    assert!(single.v1.iter().all(|component| component.is_finite()));
    assert!(single.v2.iter().all(|component| component.is_finite()));
    assert!(single.diagnostics.iters > 0);
    assert!(single.diagnostics.x.is_finite());
}

#[test]
fn invalid_request_returns_structured_error() {
    let request = LambertRequest {
        r1: [0.0, 0.0, 0.0],
        r2: [0.0, 7000.0, 0.0],
        tof: 1000.0,
        mu: MU_EARTH,
        way: TransferWayInput::Short,
        max_revs: 0,
    };

    let error = solve_lambert_request(request).unwrap_err();

    // Caller can pattern-match the variant directly — no string parsing.
    assert!(matches!(
        error,
        LambertErrorOutput::DegeneratePositionVector {
            position: PositionOutput::R1,
            ..
        }
    ));
}

#[test]
fn non_finite_input_carries_typed_parameter() {
    let request = LambertRequest {
        r1: [7000.0, f64::INFINITY, 0.0],
        r2: [0.0, 7000.0, 0.0],
        tof: 1000.0,
        mu: MU_EARTH,
        way: TransferWayInput::Short,
        max_revs: 0,
    };

    let error = solve_lambert_request(request).unwrap_err();

    assert!(matches!(
        error,
        LambertErrorOutput::NonFiniteInput {
            parameter: NonFiniteParameterOutput::R1Y,
            ..
        }
    ));
}
