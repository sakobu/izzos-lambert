//! TOF-regime scenarios: near-parabolic (Battin) dispatch and hyperbolic
//! transfers. Verified by Kepler round-trip and energy sign rather than by
//! introspecting the kernel's internal `x` (which is no longer surfaced).

use glam::DVec3;
use lambert_izzo_test_support::bodies::MU_EARTH;
use lambert_izzo_test_support::kepler::propagate as kepler_propagate;

use super::vec_sub_norm;
use crate::{LambertInput, RevolutionBudget, TransferWay, lambert};

#[test]
fn battin_regime_near_parabolic() {
    // GTO-like 90° transfer with TOF tuned so the converged x lands inside
    // the |x − 1| ≤ BATTIN_THRESHOLD (= 0.01) band, exercising the
    // hypergeometric series formulation in tof::x_to_tof_battin (Eq. 20).
    // A regression in the dispatch shows up as a Kepler round-trip blow-up.
    let mu = MU_EARTH;
    let r1 = [7000.0, 0.0, 0.0];
    let r2 = [0.0, 42_000.0, 0.0];
    let tof = 7200.0;
    let input = LambertInput {
        r1,
        r2,
        tof,
        mu,
        way: TransferWay::Short,
        revolutions: RevolutionBudget::SingleOnly,
    };
    let sols = lambert(&input).unwrap();
    assert!(sols.diagnostics.single.iters > 0);
    let r2_prop = kepler_propagate(r1, sols.single.v1, tof, mu);
    let err = vec_sub_norm(r2_prop, r2);
    assert!(err < 1e-3, "Battin round-trip err = {err} km");
}

#[test]
fn hyperbolic_transfer() {
    // Fast LEO → 200 000 km transfer; required v1 exceeds Earth escape,
    // landing in the hyperbolic branch (x > 1, positive specific energy).
    let mu = MU_EARTH;
    let r1 = [7000.0, 0.0, 0.0];
    let r2 = [0.0, 200_000.0, 0.0];
    let tof = 30_000.0;
    let input = LambertInput {
        r1,
        r2,
        tof,
        mu,
        way: TransferWay::Short,
        revolutions: RevolutionBudget::SingleOnly,
    };
    let sols = lambert(&input).unwrap();
    let v1 = DVec3::from_array(sols.single.v1);
    let energy = 0.5 * v1.dot(v1) - mu / DVec3::from_array(r1).length();
    assert!(
        energy > 0.0,
        "expected positive specific energy, got {energy}"
    );
    let r2_prop = kepler_propagate(r1, sols.single.v1, tof, mu);
    let err = vec_sub_norm(r2_prop, r2);
    assert!(err < 1e-3, "hyperbolic round-trip err = {err} km");
}
