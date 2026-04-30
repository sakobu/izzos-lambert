//! Input-validation error variants. Each test exercises one structured
//! error variant and asserts the typed payload (no string parsing).

use lambert_izzo_test_support::bodies::MU_EARTH;

use crate::{
    LambertError, LambertInput, NonFiniteParameter, Position, RevolutionBudget, TransferWay,
    constants, lambert,
};

fn valid() -> LambertInput {
    LambertInput {
        r1: [7000.0, 0.0, 0.0],
        r2: [0.0, 7000.0, 0.0],
        tof: 1000.0,
        mu: MU_EARTH,
        way: TransferWay::Short,
        revolutions: RevolutionBudget::SingleOnly,
    }
}

#[test]
fn errors_on_non_positive_tof() {
    let input = LambertInput { tof: 0.0, ..valid() };
    let err = lambert(&input).unwrap_err();
    assert!(matches!(err, LambertError::NonPositiveTimeOfFlight { tof } if tof == 0.0));
}

#[test]
fn errors_on_zero_position_vector() {
    let input = LambertInput { r1: [0.0, 0.0, 0.0], ..valid() };
    let err = lambert(&input).unwrap_err();
    assert!(matches!(
        err,
        LambertError::DegeneratePositionVector {
            position: Position::R1,
            ..
        }
    ));
}

#[test]
fn errors_on_colinear_geometry() {
    let input = LambertInput { r2: [14_000.0, 0.0, 0.0], ..valid() };
    let err = lambert(&input).unwrap_err();
    assert!(matches!(err, LambertError::CollinearGeometry { .. }));
}

#[test]
fn errors_on_transfer_angle_exactly_zero() {
    // r1 ∥ r2, same direction — sin(θ) = 0 exactly, contract is
    // |sin_angle| < COLINEARITY_TOL.
    let input = LambertInput { r2: [21_000.0, 0.0, 0.0], ..valid() };
    let err = lambert(&input).unwrap_err();
    assert!(
        matches!(
            err,
            LambertError::CollinearGeometry { sin_angle }
                if sin_angle.abs() < constants::COLINEARITY_TOL
        ),
        "expected CollinearGeometry near 0, got {err:?}",
    );
}

#[test]
fn errors_on_transfer_angle_exactly_pi() {
    // r1 anti-parallel to r2 — sin(θ) = 0 exactly, θ = π. Contract is
    // |sin_angle| < COLINEARITY_TOL.
    let input = LambertInput { r2: [-7000.0, 0.0, 0.0], ..valid() };
    let err = lambert(&input).unwrap_err();
    assert!(
        matches!(
            err,
            LambertError::CollinearGeometry { sin_angle }
                if sin_angle.abs() < constants::COLINEARITY_TOL
        ),
        "expected CollinearGeometry near π, got {err:?}",
    );
}

#[test]
fn errors_on_non_positive_mu() {
    let input = LambertInput { mu: 0.0, ..valid() };
    let err = lambert(&input).unwrap_err();
    assert!(matches!(err, LambertError::NonPositiveMu { mu } if mu == 0.0));
}

#[test]
fn errors_on_non_finite_inputs() {
    let bad_tof = LambertInput { tof: f64::NAN, ..valid() };
    let err = lambert(&bad_tof).unwrap_err();
    assert!(matches!(
        err,
        LambertError::NonFiniteInput {
            parameter: NonFiniteParameter::Tof,
            ..
        }
    ));

    let bad_r1 = LambertInput { r1: [7000.0, f64::INFINITY, 0.0], ..valid() };
    let err = lambert(&bad_r1).unwrap_err();
    assert!(matches!(
        err,
        LambertError::NonFiniteInput {
            parameter: NonFiniteParameter::R1Y,
            ..
        }
    ));
}
