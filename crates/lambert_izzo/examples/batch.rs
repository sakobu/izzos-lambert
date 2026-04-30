//! Parallel batch example. Run with:
//! `cargo run --release --example batch --features rayon`
//!
//! Demonstrates `lambert_par` over 10_000 randomized Earth-scale inputs,
//! reporting wall-clock throughput (calls/s) and the mean Householder
//! iteration count across successful single-rev solves.

use std::time::Instant;

use lambert_izzo::{RevolutionBudget, lambert_par};
use lambert_izzo_test_support::bodies::MU_EARTH;
use lambert_izzo_test_support::random_inputs::{Spec, WayStrategy, generate};
use rayon::iter::ParallelIterator;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let inputs = generate(&Spec {
        n: 10_000,
        seed: 0xA17C_0DE5,
        radius_range: (3500.0, 28_000.0),
        tof_range: (100.0, 50_000.0),
        mu: MU_EARTH,
        way: WayStrategy::Random,
        revolutions: RevolutionBudget::SingleOnly,
    });

    let n = u32::try_from(inputs.len())?;

    let start = Instant::now();
    let (successes, total_iters) = lambert_par(&inputs)
        .filter_map(Result::ok)
        .map(|sols| sols.diagnostics.single.iters)
        .fold(
            || (0u32, 0u32),
            |(c, sum), iters| (c + 1, sum + iters),
        )
        .reduce(|| (0u32, 0u32), |a, b| (a.0 + b.0, a.1 + b.1));
    let elapsed = start.elapsed();

    let throughput = f64::from(n) / elapsed.as_secs_f64();
    let mean_iters = f64::from(total_iters) / f64::from(successes);

    println!("=== lambert_par batch ({n} inputs, Earth scale) ===");
    println!("  successes: {successes} / {n}");
    println!("  wall time: {elapsed:.2?}");
    println!("  throughput: {throughput:.0} calls/s");
    println!("  mean Householder iters: {mean_iters:.3}");

    Ok(())
}
