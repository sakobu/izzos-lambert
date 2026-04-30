//! Interop tests: `nalgebra` / `glam` round-trips through the `[f64; 3]`
//! public surface, and `serde_json` round-trips when the `serde` feature
//! is on.

use core::f64::consts::PI;

use lambert_izzo_test_support::bodies::MU_EARTH;

use crate::{LambertInput, RevolutionBudget, TransferWay, lambert};

#[test]
fn interop_with_nalgebra_and_glam_round_trips() {
    let mu = MU_EARTH;
    let r1_na = nalgebra::Vector3::new(7000.0, 0.0, 0.0);
    let r2_glam = glam::DVec3::new(0.0, 7000.0, 0.0);
    let r1: [f64; 3] = r1_na.into();
    let r2: [f64; 3] = r2_glam.to_array();
    let tof = PI / 2.0 * (7000.0_f64.powi(3) / mu).sqrt();
    let input = LambertInput {
        r1,
        r2,
        tof,
        mu,
        way: TransferWay::Short,
        revolutions: RevolutionBudget::SingleOnly,
    };
    let sols = lambert(&input).unwrap();
    let v1_back_na: nalgebra::Vector3<f64> = sols.single.v1.into();
    let v2_back_glam = glam::DVec3::from_array(sols.single.v2);
    assert!(v1_back_na.iter().all(|c| c.is_finite()));
    assert!(v2_back_glam.to_array().iter().all(|c| c.is_finite()));
}

#[cfg(feature = "serde")]
#[test]
fn serde_json_round_trip_preserves_solutions_and_errors() {
    use crate::{LambertError, LambertSolutions};

    // Solutions round-trip.
    let mu = MU_EARTH;
    let r1 = [8000.0, 0.0, 0.0];
    let r2 = [5600.0, 5600.0, 0.0];
    let period = 2.0 * PI * (8000.0_f64.powi(3) / mu).sqrt();
    let input = LambertInput {
        r1,
        r2,
        tof: 5.0 * period,
        mu,
        way: TransferWay::Short,
        revolutions: RevolutionBudget::try_up_to(2).unwrap(),
    };
    let sols = lambert(&input).unwrap();
    let json = serde_json::to_string(&sols).unwrap();
    let back: LambertSolutions = serde_json::from_str(&json).unwrap();
    assert_eq!(sols, back);
    assert!(!sols.multi.is_empty(), "test should exercise multi-rev branches");

    // Error round-trip — discriminated union via the default external tag.
    let err = LambertError::CollinearGeometry { sin_angle: 1e-20 };
    let err_json = serde_json::to_string(&err).unwrap();
    assert!(err_json.contains("CollinearGeometry"));
    let err_back: LambertError = serde_json::from_str(&err_json).unwrap();
    assert_eq!(err, err_back);
}
