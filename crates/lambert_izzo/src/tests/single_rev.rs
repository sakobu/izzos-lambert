//! Single-revolution scenarios: circular LEO, heliocentric Hohmann,
//! random-input Kepler round-trip.

use core::f64::consts::PI;

use glam::DVec3;
use lambert_izzo_test_support::bodies::{AU, MU_EARTH, MU_SUN};
use lambert_izzo_test_support::kepler::propagate as kepler_propagate;
use lambert_izzo_test_support::rand_unit_vec;

use super::{approx, vec_sub_norm};
use crate::{LambertInput, RevolutionBudget, TransferWay, lambert};

#[test]
fn quarter_circle_leo() {
    // 90° transfer along a circular LEO at r = 7000 km.
    let r = 7000.0;
    let mu = MU_EARTH;
    let v_circ = (mu / r).sqrt();
    let period = 2.0 * PI * (r.powi(3) / mu).sqrt();

    let input = LambertInput {
        r1: [r, 0.0, 0.0],
        r2: [0.0, r, 0.0],
        tof: period / 4.0,
        mu,
        way: TransferWay::Short,
        revolutions: RevolutionBudget::SingleOnly,
    };
    let sols = lambert(&input).unwrap();
    assert!(sols.multi.is_empty());
    assert!(vec_sub_norm(sols.single.v1, [0.0, v_circ, 0.0]) < 1e-9);
    assert!(vec_sub_norm(sols.single.v2, [-v_circ, 0.0, 0.0]) < 1e-9);
}

#[test]
fn long_way_quarter_circle_leo() {
    // 270° transfer along the same circular LEO.
    let r = 7000.0;
    let mu = MU_EARTH;
    let v_circ = (mu / r).sqrt();
    let period = 2.0 * PI * (r.powi(3) / mu).sqrt();

    let input = LambertInput {
        r1: [r, 0.0, 0.0],
        r2: [0.0, r, 0.0],
        tof: 3.0 * period / 4.0,
        mu,
        way: TransferWay::Long,
        revolutions: RevolutionBudget::SingleOnly,
    };
    let sols = lambert(&input).unwrap();
    assert!(sols.multi.is_empty());
    assert!(vec_sub_norm(sols.single.v1, [0.0, -v_circ, 0.0]) < 1e-9);
    assert!(vec_sub_norm(sols.single.v2, [v_circ, 0.0, 0.0]) < 1e-9);
}

#[test]
fn earth_mars_hohmann() {
    // Heliocentric Hohmann transfer Earth (1 AU) → Mars (1.524 AU).
    let mu = MU_SUN;
    let r1n = AU;
    let r2n = 1.524 * AU;
    let a = f64::midpoint(r1n, r2n);
    let tof = PI * (a.powi(3) / mu).sqrt();

    let input = LambertInput {
        r1: [r1n, 0.0, 0.0],
        // 1 km off-plane to dodge the colinearity edge case.
        r2: [-r2n, 1.0, 0.0],
        tof,
        mu,
        way: TransferWay::Short,
        revolutions: RevolutionBudget::SingleOnly,
    };
    let sols = lambert(&input).unwrap();
    // Periapsis velocity: vis-viva at r = 1 AU on the transfer ellipse.
    let v_peri = (mu * (2.0 / r1n - 1.0 / a)).sqrt();
    // Tolerance ~1 m/s — the off-plane perturbation moves things slightly.
    assert!(approx(DVec3::from_array(sols.single.v1).length(), v_peri, 1e-3));
}

#[test]
fn kepler_roundtrip_random_single_rev() {
    use rand::{Rng, SeedableRng};
    use rand_chacha::ChaCha20Rng;
    use rand_distr::Uniform;

    let mu = MU_EARTH;
    let mut rng = ChaCha20Rng::seed_from_u64(0xC0FF_EE42);
    let radius = Uniform::new(3500.0, 28_000.0);
    let tof_dist = Uniform::new(100.0, 50_000.0);

    let mut max_rel_err = 0.0_f64;
    let mut good_count = 0_u32;
    let mut lambert_ok = 0_u32;
    for _ in 0..1000 {
        let (r1_mag, r2_mag) = (rng.sample(radius), rng.sample(radius));
        let r1 = (DVec3::from_array(rand_unit_vec(&mut rng)) * r1_mag).to_array();
        let r2 = (DVec3::from_array(rand_unit_vec(&mut rng)) * r2_mag).to_array();
        let tof = rng.sample(tof_dist);
        let way = if rng.gen_bool(0.5) {
            TransferWay::Long
        } else {
            TransferWay::Short
        };
        let input = LambertInput {
            r1,
            r2,
            tof,
            mu,
            way,
            revolutions: RevolutionBudget::SingleOnly,
        };
        let Ok(sols) = lambert(&input) else {
            continue;
        };
        lambert_ok += 1;
        let r2_prop = kepler_propagate(r1, sols.single.v1, tof, mu);
        let rel = vec_sub_norm(r2_prop, r2) / DVec3::from_array(r2).length();
        if rel.is_finite() {
            max_rel_err = max_rel_err.max(rel);
            good_count += 1;
        }
    }
    assert!(lambert_ok > 950, "too many Lambert failures: {lambert_ok}/1000");
    assert!(good_count > 500, "too few converged round-trips: {good_count}/{lambert_ok}");
    assert!(max_rel_err < 1e-6, "max relative round-trip err = {max_rel_err:.3e} over {good_count} trials");
}
