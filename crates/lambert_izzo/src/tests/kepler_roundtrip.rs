//! Focused single-fixture Kepler round-trip check. The random sweeps live
//! under `single_rev` / `multi_rev`; this test pins one specific geometry
//! so a regression localizes immediately.

use lambert_izzo_test_support::bodies::MU_EARTH;
use lambert_izzo_test_support::kepler::propagate as kepler_propagate;

use super::vec_sub_norm;
use crate::{LambertInput, RevolutionBudget, TransferWay, lambert};

#[test]
fn round_trip_kepler_check_single_rev() {
    // Propagate v1 with a universal-variable Kepler integrator and confirm
    // we land within numerical tolerance of r2.
    let mu = MU_EARTH;
    let r1 = [10_500.0, 1400.0, 700.0];
    let r2 = [-2800.0, 9100.0, -1400.0];
    let tof = 4500.0;
    let input = LambertInput {
        r1,
        r2,
        tof,
        mu,
        way: TransferWay::Short,
        revolutions: RevolutionBudget::SingleOnly,
    };
    let sols = lambert(&input).unwrap();
    let v1 = sols.single.v1;
    let r2_prop = kepler_propagate(r1, v1, tof, mu);
    let err = vec_sub_norm(r2_prop, r2);
    assert!(err < 1e-6, "kepler-roundtrip err = {err} km");
}
