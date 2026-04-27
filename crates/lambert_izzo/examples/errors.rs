//! Worked tour of `LambertError` — when each variant fires and how to
//! recover. Run with `cargo run --release --example errors`.

use core::f64::consts::PI;

use lambert_izzo::{
    LambertError, NonFiniteParameter, RevolutionBudget, TransferWay, lambert,
    solve_with_diagnostics,
};

const MU_EARTH_KM3_S2: f64 = 398_600.441_8;

fn main() {
    println!("=== 1. Colinear geometry — perturb to recover ===");
    // r1 and r2 both on the +X axis: transfer plane undefined.
    let r1_km = [7000.0, 0.0, 0.0];
    let r2_km = [14_000.0, 0.0, 0.0];
    match lambert(
        r1_km,
        r2_km,
        1500.0,
        MU_EARTH_KM3_S2,
        TransferWay::Short,
        RevolutionBudget::SingleOnly,
    ) {
        Err(LambertError::CollinearGeometry { sin_angle }) => {
            println!(
                "  rejected: sin(θ) = {sin_angle:.3e} (below COLINEARITY_TOL).\n  \
                 Recovery: perturb r2 by 1 km off-plane."
            );
            // Re-solve with a tiny off-plane component.
            let r2_perturbed = [r2_km[0], 1.0, 0.0];
            let sols = lambert(
                r1_km,
                r2_perturbed,
                1500.0,
                MU_EARTH_KM3_S2,
                TransferWay::Short,
                RevolutionBudget::SingleOnly,
            )
            .expect("perturbed geometry should converge");
            println!(
                "  perturbed v1 = [{:+.4}, {:+.4}, {:+.4}] km/s",
                sols.single.v1_km_s[0], sols.single.v1_km_s[1], sols.single.v1_km_s[2],
            );
        }
        Err(e) => println!("  unexpected error: {e}"),
        Ok(_) => println!("  unexpected success"),
    }

    println!("\n=== 2. Non-finite input — parameter is typed, not stringly ===");
    let bad_r1 = [7000.0, f64::INFINITY, 0.0];
    let r2_km = [0.0, 7000.0, 0.0];
    match lambert(
        bad_r1,
        r2_km,
        1500.0,
        MU_EARTH_KM3_S2,
        TransferWay::Short,
        RevolutionBudget::SingleOnly,
    ) {
        Err(LambertError::NonFiniteInput { parameter, value }) => {
            // Pattern-match on the typed enum, not a string.
            let component = match parameter {
                NonFiniteParameter::R1KmY => "r1.y",
                other => other.as_str(),
            };
            println!("  rejected: {component} = {value} (must be finite).");
        }
        other => println!("  unexpected: {other:?}"),
    }

    println!("\n=== 3. Multi-rev infeasibility — silent skip ===");
    // Earth-orbit phasing with up_to(10), but tof only large enough for M=1
    // and M=2. Higher M get dropped from the returned `multi` list.
    let r1_km = [8000.0, 0.0, 0.0];
    let r2_km = [5600.0, 5600.0, 0.0];
    let period_s = 2.0 * PI * (8000.0_f64.powi(3) / MU_EARTH_KM3_S2).sqrt();
    let tof_s = 3.0 * period_s; // budget allows M up to ~3, but T_min cuts it short
    let sols = lambert(
        r1_km,
        r2_km,
        tof_s,
        MU_EARTH_KM3_S2,
        TransferWay::Short,
        RevolutionBudget::up_to(10),
    )
    .expect("phasing should converge");
    println!(
        "  asked for up to M=10, got {} multi-rev pair(s) — solver dropped \
         branches whose T_min exceeded tof.",
        sols.multi.len()
    );
    for pair in &sols.multi {
        println!(
            "    M={}: long-period |v1|={:.3} km/s, short-period |v1|={:.3} km/s",
            pair.n_revs,
            (pair.long_period.v1_km_s[0].powi(2)
                + pair.long_period.v1_km_s[1].powi(2)
                + pair.long_period.v1_km_s[2].powi(2))
            .sqrt(),
            (pair.short_period.v1_km_s[0].powi(2)
                + pair.short_period.v1_km_s[1].powi(2)
                + pair.short_period.v1_km_s[2].powi(2))
            .sqrt(),
        );
    }

    println!("\n=== 4. Near-parabolic — Battin dispatch is automatic ===");
    // GTO-like 90° transfer with a TOF tuned to land x in the Battin band.
    let r1_km = [7000.0, 0.0, 0.0];
    let r2_km = [0.0, 42_000.0, 0.0];
    let tof_s = 7200.0;
    let (sols, diag) = solve_with_diagnostics(
        r1_km,
        r2_km,
        tof_s,
        MU_EARTH_KM3_S2,
        TransferWay::Short,
        RevolutionBudget::SingleOnly,
    )
    .expect("near-parabolic should converge via Battin");
    let x = diag.single.lancaster_blanchard_x;
    println!(
        "  converged x = {x:.6}, |x−1| = {:.3e}; Battin threshold = {:.3e}.",
        (x - 1.0).abs(),
        lambert_izzo::constants::BATTIN_THRESHOLD,
    );
    println!(
        "  iters = {} (single-rev paper avg ≈ 2.1).",
        diag.single.iters,
    );
    println!(
        "  trajectory v1 = [{:+.4}, {:+.4}, {:+.4}] km/s",
        sols.single.v1_km_s[0], sols.single.v1_km_s[1], sols.single.v1_km_s[2],
    );

    println!("\n=== 5. Non-positive scalar — distinct error variants ===");
    for (label, err) in [
        (
            "tof_s = 0",
            lambert(
                r1_km,
                r2_km,
                0.0,
                MU_EARTH_KM3_S2,
                TransferWay::Short,
                RevolutionBudget::SingleOnly,
            ),
        ),
        (
            "mu_km3_s2 = -1",
            lambert(
                r1_km,
                r2_km,
                1500.0,
                -1.0,
                TransferWay::Short,
                RevolutionBudget::SingleOnly,
            ),
        ),
    ] {
        match err {
            Err(LambertError::NonPositiveTimeOfFlight { tof_s }) => {
                println!("  {label}: NonPositiveTimeOfFlight(tof_s = {tof_s})");
            }
            Err(LambertError::NonPositiveMu { mu_km3_s2 }) => {
                println!("  {label}: NonPositiveMu(mu_km3_s2 = {mu_km3_s2})");
            }
            Err(other) => println!("  {label}: unexpected variant {other:?}"),
            Ok(_) => println!("  {label}: unexpected success"),
        }
    }
}
