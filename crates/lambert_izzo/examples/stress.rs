//! Stress test: re-create the experiment from §5 of Izzo (2014) at Earth scale.
//!
//! Verifies the returned `(v1, v2)` satisfy orbital conservation laws —
//! vis-viva energy and angular momentum — which doesn't depend on any
//! secondary numerical integrator. Clean signal on the solver itself.
//!
//! Run with `cargo run --release --example stress`.

use glam::DVec3;
use lambert_izzo::{
    LambertSolution, RevolutionBudget, TransferWay, lambert, solve_with_diagnostics,
};
use lambert_izzo_test_support::bodies::MU_EARTH;
use lambert_izzo_test_support::rand_unit_vec;
use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha20Rng;
use rand_distr::Uniform;

fn check_conservation(r1: [f64; 3], s: LambertSolution, r2: [f64; 3], mu: f64) -> (f64, f64) {
    let r1v = DVec3::from_array(r1);
    let r2v = DVec3::from_array(r2);
    let v1 = DVec3::from_array(s.v1);
    let v2 = DVec3::from_array(s.v2);
    let e1 = 0.5 * v1.dot(v1) - mu / r1v.length();
    let e2 = 0.5 * v2.dot(v2) - mu / r2v.length();
    let e_avg = 0.5 * (e1.abs() + e2.abs()).max(1e-12);
    let e_rel = (e1 - e2).abs() / e_avg;
    let h1 = r1v.cross(v1);
    let h2 = r2v.cross(v2);
    let h_avg = 0.5 * (h1.length() + h2.length()).max(1e-12);
    let h_rel = (h1 - h2).length() / h_avg;
    (e_rel, h_rel)
}

fn scaled_unit(rng: &mut ChaCha20Rng, scalar: f64) -> [f64; 3] {
    let u = rand_unit_vec(rng);
    [u[0] * scalar, u[1] * scalar, u[2] * scalar]
}

fn main() {
    let n_trials = 100_000_u32;
    let mu = MU_EARTH;

    // Single-revolution sweep: r ∈ [3500, 28000] km, tof ∈ [100, 50000] s.
    {
        let mut rng = ChaCha20Rng::seed_from_u64(0x00C0_FFEE);
        let radius = Uniform::new(3500.0, 28_000.0);
        let tof_dist = Uniform::new(100.0, 50_000.0);
        let mut iters_hist = [0_u32; 16];
        let mut total_iters = 0_u32;
        let mut count = 0_u32;
        let mut nonconvergence = 0_u32;
        let mut max_e_rel = 0.0_f64;
        let mut max_h_rel = 0.0_f64;
        for _ in 0..n_trials {
            let (r1_mag, r2_mag) = (rng.sample(radius), rng.sample(radius));
            let r1 = scaled_unit(&mut rng, r1_mag);
            let r2 = scaled_unit(&mut rng, r2_mag);
            let tof = rng.sample(tof_dist);
            let way = if rng.gen_bool(0.5) {
                TransferWay::Long
            } else {
                TransferWay::Short
            };
            match solve_with_diagnostics(r1, r2, tof, mu, way, RevolutionBudget::SingleOnly) {
                Ok((sols, diag)) => {
                    iters_hist[diag.single.iters as usize] += 1; // u32 → usize: always safe
                    total_iters += diag.single.iters;
                    let (e_rel, h_rel) = check_conservation(r1, sols.single, r2, mu);
                    if e_rel.is_finite() {
                        max_e_rel = max_e_rel.max(e_rel);
                    }
                    if h_rel.is_finite() {
                        max_h_rel = max_h_rel.max(h_rel);
                    }
                    count += 1;
                }
                Err(_) => nonconvergence += 1,
            }
        }
        println!("=== Single revolution ({n_trials} trials, Earth scale) ===");
        println!("  successful: {count}, non-convergence: {nonconvergence}");
        println!(
            "  avg iterations: {:.3}  (paper: 2.1)",
            f64::from(total_iters) / f64::from(count)
        );
        print!("  iter histogram: ");
        for (i, c) in iters_hist.iter().enumerate().filter(|&(_, c)| *c > 0) {
            print!("{i}={c} ");
        }
        println!();
        println!("  max relative energy mismatch: {max_e_rel:.2e}");
        println!("  max relative ang.mom. mismatch: {max_h_rel:.2e}");
    }

    // Multi-revolution sweep: r ∈ [5600, 10500] km, tof ∈ [10000, 250000] s.
    {
        let mut rng = ChaCha20Rng::seed_from_u64(0xDEAD_BEEF);
        let radius = Uniform::new(5600.0, 10_500.0);
        let tof_dist = Uniform::new(10_000.0, 250_000.0);
        let mut iters_hist = [0_u32; 16];
        let mut total_iters = 0_u32;
        let mut total_branches = 0_u32;
        let mut nonconvergence = 0_u32;
        let mut max_e_rel = 0.0_f64;
        let mut max_h_rel = 0.0_f64;
        for _ in 0..n_trials {
            let (r1_mag, r2_mag) = (rng.sample(radius), rng.sample(radius));
            let r1 = scaled_unit(&mut rng, r1_mag);
            let r2 = scaled_unit(&mut rng, r2_mag);
            let tof = rng.sample(tof_dist);
            match solve_with_diagnostics(
                r1,
                r2,
                tof,
                mu,
                TransferWay::Short,
                RevolutionBudget::up_to(5),
            ) {
                Ok((sols, diag)) => {
                    for (pair, dpair) in sols.multi.iter().zip(diag.multi.iter()) {
                        for (branch, dbranch) in [
                            (pair.long_period, dpair.long_period),
                            (pair.short_period, dpair.short_period),
                        ] {
                            iters_hist[dbranch.iters as usize] += 1; // u32 → usize: always safe
                            total_iters += dbranch.iters;
                            let (e_rel, h_rel) = check_conservation(r1, branch, r2, mu);
                            if e_rel.is_finite() {
                                max_e_rel = max_e_rel.max(e_rel);
                            }
                            if h_rel.is_finite() {
                                max_h_rel = max_h_rel.max(h_rel);
                            }
                            total_branches += 1;
                        }
                    }
                }
                Err(_) => nonconvergence += 1,
            }
        }
        println!("\n=== Multi revolution ({n_trials} trials, M up to 5, Earth scale) ===");
        println!("  multi-rev branches: {total_branches}, non-convergence: {nonconvergence}");
        if total_branches > 0 {
            println!(
                "  avg iterations: {:.3}  (paper: 3.3)",
                f64::from(total_iters) / f64::from(total_branches)
            );
            print!("  iter histogram: ");
            for (i, c) in iters_hist.iter().enumerate().filter(|&(_, c)| *c > 0) {
                print!("{i}={c} ");
            }
            println!();
            println!("  max relative energy mismatch: {max_e_rel:.2e}");
            println!("  max relative ang.mom. mismatch: {max_h_rel:.2e}");
        }
    }

    // Sanity-check the bare lambert() path.
    let r1 = [7000.0, 0.0, 0.0];
    let r2 = [0.0, 7000.0, 0.0];
    let tof = core::f64::consts::PI / 2.0 * (7000.0_f64.powi(3) / mu).sqrt();
    let _ = lambert(
        r1,
        r2,
        tof,
        mu,
        TransferWay::Short,
        RevolutionBudget::SingleOnly,
    )
    .expect("trivial LEO Lambert call should succeed");
}
