//! Smoke test for the `test-utils` feature: round-trip a Lambert solution
//! through the publicly-exported `kepler_propagate` and verify it lands
//! back at `r2`. Compiled only when `--features test-utils` is enabled.

#![cfg(feature = "test-utils")]

use lambert_izzo::{RevolutionBudget, TransferWay, lambert, test_utils::kepler_propagate};

const MU_EARTH: f64 = 398_600.441_8;

fn vec_sub_norm(a: [f64; 3], b: [f64; 3]) -> f64 {
    let d = [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
    (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
}

#[test]
fn kepler_propagate_round_trips_a_lambert_solution() {
    let r1 = [10_500.0, 1400.0, 700.0];
    let r2 = [-2800.0, 9100.0, -1400.0];
    let tof = 4500.0;

    let sols = lambert(
        r1,
        r2,
        tof,
        MU_EARTH,
        TransferWay::Short,
        RevolutionBudget::SingleOnly,
    )
    .expect("Lambert should converge");

    let r2_prop = kepler_propagate(r1, sols.single.v1, tof, MU_EARTH);
    let err = vec_sub_norm(r2_prop, r2);
    assert!(err < 1e-6, "kepler-roundtrip err = {err} km");
}
