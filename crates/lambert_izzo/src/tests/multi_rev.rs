//! Multi-revolution scenarios: branching, ordering, distinct branches,
//! random-input Kepler round-trip.

use core::f64::consts::PI;

use glam::DVec3;
use lambert_izzo_test_support::bodies::MU_EARTH;
use lambert_izzo_test_support::kepler::propagate as kepler_propagate;
use lambert_izzo_test_support::rand_unit_vec;

use super::vec_sub_norm;
use crate::{LambertInput, LambertSolution, RevolutionBudget, TransferWay, lambert};

/// Long-tof Earth-orbit phasing — admits multiple multi-rev branches.
/// Shared across the three scenarios that exercise branching, ordering,
/// and round-trip distinctness on the same geometry.
fn phasing_input() -> LambertInput {
    let r = 8000.0_f64;
    let period = 2.0 * PI * (r.powi(3) / MU_EARTH).sqrt();
    LambertInput {
        r1: [r, 0.0, 0.0],
        r2: [5600.0, 5600.0, 0.0],
        tof: 5.0 * period,
        mu: MU_EARTH,
        way: TransferWay::Short,
        revolutions: RevolutionBudget::try_up_to(3).unwrap(),
    }
}

#[test]
fn multi_rev_branches() {
    let input = phasing_input();
    let sols = lambert(&input).unwrap();
    assert!(!sols.multi.is_empty(), "no multi-rev pairs returned");
    let r1n = DVec3::from_array(input.r1).length();
    for pair in &sols.multi {
        for s in [pair.long_period, pair.short_period] {
            let v1 = DVec3::from_array(s.v1);
            let energy = 0.5 * v1.dot(v1) - input.mu / r1n;
            assert!(energy.is_finite());
        }
    }
}

#[test]
fn multi_rev_branches_distinct() {
    // Long-period vs short-period branches must produce distinct trajectories
    // and both must round-trip via the universal-variable Kepler propagator.
    let input = phasing_input();
    let sols = lambert(&input).unwrap();

    assert!(!sols.multi.is_empty(), "no multi-rev pairs found");
    assert_eq!(sols.multi.len(), sols.diagnostics.multi.len());

    for (pair, diag_pair) in sols.multi.iter().zip(sols.diagnostics.multi.iter()) {
        assert_eq!(pair.n_revs, diag_pair.n_revs);
        // Long-period and short-period must produce different trajectories.
        assert!(
            vec_sub_norm(pair.long_period.v1, pair.short_period.v1) > 1e-6,
            "M = {}: long-period and short-period v1 are identical",
            pair.n_revs,
        );
        for branch in [pair.long_period, pair.short_period] {
            let r2_prop = kepler_propagate(input.r1, branch.v1, input.tof, input.mu);
            let err = vec_sub_norm(r2_prop, input.r2);
            assert!(err < 1e-3, "M = {} branch round-trip err = {err} km", pair.n_revs);
        }
    }
}

#[test]
fn solution_ordering_contract() {
    let input = phasing_input();
    let sols = lambert(&input).unwrap();

    let mut prev_m = 0_u32;
    for pair in &sols.multi {
        assert!(pair.n_revs > prev_m, "M strictly ascending across pairs");
        prev_m = pair.n_revs;
    }
}

#[test]
fn kepler_roundtrip_random_multi_rev() {
    use rand::{Rng, SeedableRng};
    use rand_chacha::ChaCha20Rng;
    use rand_distr::Uniform;

    let mu = MU_EARTH;
    let mut rng = ChaCha20Rng::seed_from_u64(0xBEEF_DEAD);
    let radius = Uniform::new(5600.0, 10_500.0);
    let tof_dist = Uniform::new(10_000.0, 250_000.0);

    let mut max_rel_err = 0.0_f64;
    let mut branches = 0_u32;
    let mut good_count = 0_u32;
    for _ in 0..500 {
        let (r1_mag, r2_mag) = (rng.sample(radius), rng.sample(radius));
        let r1 = (DVec3::from_array(rand_unit_vec(&mut rng)) * r1_mag).to_array();
        let r2 = (DVec3::from_array(rand_unit_vec(&mut rng)) * r2_mag).to_array();
        let tof = rng.sample(tof_dist);
        let input = LambertInput {
            r1,
            r2,
            tof,
            mu,
            way: TransferWay::Short,
            revolutions: RevolutionBudget::try_up_to(3).unwrap(),
        };
        let Ok(sols) = lambert(&input) else {
            continue;
        };
        let mut iter_branch = |s: LambertSolution| {
            let r2_prop = kepler_propagate(r1, s.v1, tof, mu);
            let rel = vec_sub_norm(r2_prop, r2) / DVec3::from_array(r2).length();
            if rel.is_finite() {
                max_rel_err = max_rel_err.max(rel);
                good_count += 1;
            }
            branches += 1;
        };
        iter_branch(sols.single);
        for pair in &sols.multi {
            iter_branch(pair.long_period);
            iter_branch(pair.short_period);
        }
    }
    assert!(branches > 500, "expected branches > 500, got {branches}");
    assert!(good_count > 300, "too few converged round-trips: {good_count}/{branches}");
    assert!(max_rel_err < 1e-5, "max relative round-trip err = {max_rel_err:.3e} over {good_count} branches");
}

#[test]
fn max_feasible_revs_reflects_silent_skip() {
    // Earth-orbit phasing with up_to(10), but tof only large enough for a
    // few multi-rev branches. `max_feasible_revs` must report the highest
    // M that actually produced a pair, below the requested cap.
    let r = 8000.0_f64;
    let period = 2.0 * PI * (r.powi(3) / MU_EARTH).sqrt();
    let input = LambertInput {
        r1: [r, 0.0, 0.0],
        r2: [5600.0, 5600.0, 0.0],
        tof: 3.0 * period,
        mu: MU_EARTH,
        way: TransferWay::Short,
        revolutions: RevolutionBudget::try_up_to(10).unwrap(),
    };
    let sols = lambert(&input).unwrap();

    assert!(!sols.multi.is_empty(), "no multi-rev pairs returned");
    assert!(
        sols.max_feasible_revs() < 10,
        "expected silent-skip below requested cap; got max_feasible_revs = {}",
        sols.max_feasible_revs(),
    );
}

#[test]
fn max_feasible_revs_zero_when_single_only() {
    // Same phasing geometry, but the budget rejects every multi-rev branch.
    let mut input = phasing_input();
    input.revolutions = RevolutionBudget::SingleOnly;
    let sols = lambert(&input).unwrap();

    assert!(sols.multi.is_empty());
    assert_eq!(sols.max_feasible_revs(), 0);
}

#[test]
fn max_feasible_revs_zero_when_no_multi_branch_feasible() {
    // Quarter-circle Earth orbit — TOF is one quarter period, far below
    // T_min for every M ≥ 1, so no multi-rev branch is feasible even
    // though the budget asked for some.
    let r = 7000.0_f64;
    let tof = (PI / 2.0) * (r.powi(3) / MU_EARTH).sqrt();
    let input = LambertInput {
        r1: [r, 0.0, 0.0],
        r2: [0.0, r, 0.0],
        tof,
        mu: MU_EARTH,
        way: TransferWay::Short,
        revolutions: RevolutionBudget::try_up_to(3).unwrap(),
    };
    let sols = lambert(&input).unwrap();

    assert!(sols.multi.is_empty());
    assert_eq!(sols.max_feasible_revs(), 0);
}
