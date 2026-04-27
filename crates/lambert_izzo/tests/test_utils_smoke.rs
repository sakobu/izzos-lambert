//! Smoke test for the `test-utils` feature: round-trip a Lambert solution
//! through the publicly-exported `kepler_propagate` and verify it lands
//! back at `r2`. Compiled only when `--features test-utils` is enabled.

#![cfg(feature = "test-utils")]

use lambert_izzo::{RevolutionBudget, TransferWay, lambert, test_utils::kepler_propagate};

const MU_EARTH_KM3_S2: f64 = 398_600.441_8;

fn vec_sub_norm(a: [f64; 3], b: [f64; 3]) -> f64 {
    let d = [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
    (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
}

#[test]
fn kepler_propagate_round_trips_a_lambert_solution() {
    let r1_km = [10_500.0, 1400.0, 700.0];
    let r2_km = [-2800.0, 9100.0, -1400.0];
    let tof_s = 4500.0;

    let sols = lambert(
        r1_km,
        r2_km,
        tof_s,
        MU_EARTH_KM3_S2,
        TransferWay::Short,
        RevolutionBudget::SingleOnly,
    )
    .expect("Lambert should converge");

    let r2_prop = kepler_propagate(r1_km, sols.single.v1_km_s, tof_s, MU_EARTH_KM3_S2);
    let err_km = vec_sub_norm(r2_prop, r2_km);
    assert!(err_km < 1e-6, "kepler-roundtrip err = {err_km} km");
}
