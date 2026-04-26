//! Stress test: re-create the experiment from §5 of Izzo (2014) at Earth scale.
//!
//! Verifies the returned `(v1_km_s, v2_km_s)` satisfy orbital conservation
//! laws — vis-viva energy and angular momentum — which doesn't depend on any
//! secondary numerical integrator. Clean signal on the solver itself.
//!
//! Run with `cargo run --release --example stress`.

use lambert_izzo::{lambert, RevolutionBudget, TransferWay};
use nalgebra::Vector3;
use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha20Rng;
use rand_distr::Uniform;

/// Earth's gravitational parameter (km³/s²).
const MU_EARTH_KM3_S2: f64 = 398_600.441_8;

fn rand_unit_vec(rng: &mut ChaCha20Rng) -> Vector3<f64> {
    let axis: Uniform<f64> = Uniform::new(-1.0, 1.0);
    loop {
        let v = Vector3::new(rng.sample(axis), rng.sample(axis), rng.sample(axis));
        let n2 = v.norm_squared();
        if n2 > 0.01 && n2 < 1.0 {
            return v / n2.sqrt();
        }
    }
}

fn check_conservation(
    r1_km: Vector3<f64>,
    v1_km_s: Vector3<f64>,
    r2_km: Vector3<f64>,
    v2_km_s: Vector3<f64>,
    mu_km3_s2: f64,
) -> (f64, f64) {
    let r1n = r1_km.norm();
    let r2n = r2_km.norm();
    let e1 = 0.5 * v1_km_s.dot(&v1_km_s) - mu_km3_s2 / r1n;
    let e2 = 0.5 * v2_km_s.dot(&v2_km_s) - mu_km3_s2 / r2n;
    let e_avg = 0.5 * (e1.abs() + e2.abs()).max(1e-12);
    let e_rel = (e1 - e2).abs() / e_avg;
    let h1 = r1_km.cross(&v1_km_s);
    let h2 = r2_km.cross(&v2_km_s);
    let h_diff = (h1 - h2).norm();
    let h_avg = 0.5 * (h1.norm() + h2.norm()).max(1e-12);
    let h_rel = h_diff / h_avg;
    (e_rel, h_rel)
}

fn main() {
    let n_trials = 100_000_u32;
    let mu = MU_EARTH_KM3_S2;

    // Single-revolution sweep: r ∈ [3500, 28000] km, tof ∈ [100, 50000] s.
    {
        let mut rng = ChaCha20Rng::seed_from_u64(0x00C0_FFEE);
        let radius = Uniform::new(3500.0, 28_000.0);
        let tof = Uniform::new(100.0, 50_000.0);
        let mut iters_hist = [0_u32; 16];
        let mut total_iters = 0_u32;
        let mut count = 0_u32;
        let mut nonconvergence = 0_u32;
        let mut max_e_rel = 0.0_f64;
        let mut max_h_rel = 0.0_f64;
        for _ in 0..n_trials {
            let r1_km = rand_unit_vec(&mut rng) * rng.sample(radius);
            let r2_km = rand_unit_vec(&mut rng) * rng.sample(radius);
            let tof_s = rng.sample(tof);
            let way = if rng.gen_bool(0.5) { TransferWay::Long } else { TransferWay::Short };
            match lambert(r1_km, r2_km, tof_s, mu, way, RevolutionBudget::SingleOnly) {
                Ok(sols) => {
                    let s = sols[0];
                    iters_hist[s.iters as usize] += 1; // u32 → usize: always safe (usize ≥ 32 bits)
                    total_iters += s.iters;
                    let (e_rel, h_rel) =
                        check_conservation(r1_km, s.v1_km_s, r2_km, s.v2_km_s, mu);
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
        let tof = Uniform::new(10_000.0, 250_000.0);
        let mut iters_hist = [0_u32; 16];
        let mut total_iters = 0_u32;
        let mut total_branches = 0_u32;
        let mut nonconvergence = 0_u32;
        let mut max_e_rel = 0.0_f64;
        let mut max_h_rel = 0.0_f64;
        for _ in 0..n_trials {
            let r1_km = rand_unit_vec(&mut rng) * rng.sample(radius);
            let r2_km = rand_unit_vec(&mut rng) * rng.sample(radius);
            let tof_s = rng.sample(tof);
            match lambert(r1_km, r2_km, tof_s, mu, TransferWay::Short, RevolutionBudget::up_to(5)) {
                Ok(sols) => {
                    for s in sols.iter().filter(|s| s.n_revs > 0) {
                        iters_hist[s.iters as usize] += 1; // u32 → usize: always safe (usize ≥ 32 bits)
                        total_iters += s.iters;
                        let (e_rel, h_rel) =
                            check_conservation(r1_km, s.v1_km_s, r2_km, s.v2_km_s, mu);
                        if e_rel.is_finite() {
                            max_e_rel = max_e_rel.max(e_rel);
                        }
                        if h_rel.is_finite() {
                            max_h_rel = max_h_rel.max(h_rel);
                        }
                        total_branches += 1;
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
}
