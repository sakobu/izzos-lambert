use lambert_izzo_test_support::bodies::MU_EARTH;
use lambert_izzo_wasm::{
    BatchResult, LambertErrorOutput, LambertRequest, NonFiniteParameterOutput, PositionOutput,
    TransferWayInput, solve_lambert_batch, solve_lambert_request,
};

fn earth_quarter_circle_request() -> LambertRequest {
    LambertRequest {
        r1: [7000.0, 0.0, 0.0],
        r2: [0.0, 7000.0, 0.0],
        tof: core::f64::consts::PI / 2.0 * (7000.0_f64.powi(3) / MU_EARTH).sqrt(),
        mu: MU_EARTH,
        way: TransferWayInput::Short,
        max_revs: None,
    }
}

#[test]
fn single_rev_request_returns_js_friendly_vectors() {
    let response = solve_lambert_request(earth_quarter_circle_request()).unwrap();

    assert!(response.multi.is_empty());
    let single = response.single;
    assert!(single.v1.iter().all(|component| component.is_finite()));
    assert!(single.v2.iter().all(|component| component.is_finite()));
    assert!(response.diagnostics.single.iters > 0);
    assert!(response.diagnostics.multi.is_empty());
}

#[test]
fn invalid_request_returns_structured_error() {
    let request = LambertRequest {
        r1: [0.0, 0.0, 0.0],
        r2: [0.0, 7000.0, 0.0],
        tof: 1000.0,
        mu: MU_EARTH,
        way: TransferWayInput::Short,
        max_revs: None,
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
        max_revs: None,
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

#[test]
fn out_of_range_max_revs_returns_structured_error() {
    // 0 is below the type-enforced minimum (must be >= 1); the request fails
    // before any solver work happens.
    let request = LambertRequest {
        max_revs: Some(0),
        ..earth_quarter_circle_request()
    };
    let error = solve_lambert_request(request).unwrap_err();
    assert!(matches!(
        error,
        LambertErrorOutput::RevsOutOfRange {
            requested: 0,
            max: 32,
        }
    ));

    // Above-cap value is also rejected at request time.
    let request = LambertRequest {
        max_revs: Some(100),
        ..earth_quarter_circle_request()
    };
    let error = solve_lambert_request(request).unwrap_err();
    assert!(matches!(
        error,
        LambertErrorOutput::RevsOutOfRange {
            requested: 100,
            max: 32,
        }
    ));
}

#[test]
fn batch_request_returns_per_input_results() {
    // Batch of three: valid, degenerate r1, non-finite r1.y.
    let valid = earth_quarter_circle_request();
    let degenerate = LambertRequest {
        r1: [0.0, 0.0, 0.0],
        ..valid
    };
    let non_finite = LambertRequest {
        r1: [7000.0, f64::INFINITY, 0.0],
        ..valid
    };

    let results = solve_lambert_batch(vec![valid, degenerate, non_finite]);

    assert_eq!(results.len(), 3);
    assert!(matches!(&results[0], BatchResult::Ok { .. }));
    assert!(matches!(
        &results[1],
        BatchResult::Err {
            error: LambertErrorOutput::DegeneratePositionVector {
                position: PositionOutput::R1,
                ..
            }
        }
    ));
    assert!(matches!(
        &results[2],
        BatchResult::Err {
            error: LambertErrorOutput::NonFiniteInput {
                parameter: NonFiniteParameterOutput::R1Y,
                ..
            }
        }
    ));
}
