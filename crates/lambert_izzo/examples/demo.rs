//! Worked examples for the [`lambert_izzo`] crate.
//!
//! Run with `cargo run --release --example demo`.

use core::f64::consts::PI;
use lambert_izzo::{
    LambertSolution, RevolutionBudget, TransferWay, lambert, solve_with_diagnostics,
};

const MU_EARTH_KM3_S2: f64 = 398_600.441_8;
const MU_SUN_KM3_S2: f64 = 1.327_124_400_18e11;
const AU_KM: f64 = 1.495_978_707e8;

fn norm(v: [f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

fn print_trajectory(label: &str, s: LambertSolution) {
    println!(
        "  {label}: v1 = [{:+.6}, {:+.6}, {:+.6}] km/s  |v1| = {:.4} km/s  |v2| = {:.4} km/s",
        s.v1_km_s[0],
        s.v1_km_s[1],
        s.v1_km_s[2],
        norm(s.v1_km_s),
        norm(s.v2_km_s),
    );
}

fn main() {
    // 1. LEO → MEO Hohmann transfer (Earth-centered).
    println!("=== LEO (7000 km) → MEO (12000 km) Hohmann, short way ===");
    let r1_km = [7000.0, 0.0, 0.0];
    let r2_km = [-12_000.0, 1.0, 0.0]; // 1 km off-axis for non-collinearity
    let a_km = f64::midpoint(7000.0, 12_000.0);
    let tof_s = PI * (a_km.powi(3) / MU_EARTH_KM3_S2).sqrt();
    let sols = lambert(
        r1_km,
        r2_km,
        tof_s,
        MU_EARTH_KM3_S2,
        TransferWay::Short,
        RevolutionBudget::SingleOnly,
    )
    .expect("LEO Hohmann should converge");
    print_trajectory("single-rev", sols.single);

    // 2. Same geometry, long way around.
    println!("\n=== Same geometry, long way ===");
    let sols = lambert(
        r1_km,
        r2_km,
        tof_s,
        MU_EARTH_KM3_S2,
        TransferWay::Long,
        RevolutionBudget::SingleOnly,
    )
    .expect("long-way Hohmann should converge");
    print_trajectory("single-rev", sols.single);

    // 3. Earth → Mars heliocentric Hohmann (Sun-centered, large scale).
    println!("\n=== Earth (1 AU) → Mars (1.524 AU) heliocentric, short way ===");
    let r1_km = [AU_KM, 0.0, 0.0];
    let r2_km = [-1.524 * AU_KM, 1.0, 0.0];
    let a_km = f64::midpoint(AU_KM, 1.524 * AU_KM);
    let tof_s = PI * (a_km.powi(3) / MU_SUN_KM3_S2).sqrt();
    let sols = lambert(
        r1_km,
        r2_km,
        tof_s,
        MU_SUN_KM3_S2,
        TransferWay::Short,
        RevolutionBudget::SingleOnly,
    )
    .expect("Earth-Mars Hohmann should converge");
    print_trajectory(
        &format!("single-rev (tof = {:.2} days)", tof_s / 86400.0),
        sols.single,
    );

    // 4. Earth-orbit multi-rev phasing — show the new pair structure.
    println!("\n=== Earth-orbit phasing (long tof, M up to 3) ===");
    let r1_km = [8000.0, 0.0, 0.0];
    let r2_km = [5600.0, 5600.0, 0.0];
    let period_s = 2.0 * PI * (8000.0_f64.powi(3) / MU_EARTH_KM3_S2).sqrt();
    let tof_s = 5.0 * period_s;
    let (sols, diag) = solve_with_diagnostics(
        r1_km,
        r2_km,
        tof_s,
        MU_EARTH_KM3_S2,
        TransferWay::Short,
        RevolutionBudget::up_to(3),
    )
    .expect("multi-rev phasing should converge");
    println!(
        "  single-rev: iters={} x={:+.6}  |v1|={:.4} km/s",
        diag.single.iters,
        diag.single.lancaster_blanchard_x,
        norm(sols.single.v1_km_s),
    );
    for (pair, dpair) in sols.multi.iter().zip(diag.multi.iter()) {
        println!(
            "  M={} long-period:  iters={} x={:+.6}  |v1|={:.4} km/s",
            pair.n_revs,
            dpair.long_period.iters,
            dpair.long_period.lancaster_blanchard_x,
            norm(pair.long_period.v1_km_s),
        );
        println!(
            "  M={} short-period: iters={} x={:+.6}  |v1|={:.4} km/s",
            pair.n_revs,
            dpair.short_period.iters,
            dpair.short_period.lancaster_blanchard_x,
            norm(pair.short_period.v1_km_s),
        );
    }
}
