//! Worked examples for the [`lambert_izzo`] crate.
//!
//! Run with `cargo run --release --example demo`.

use core::f64::consts::PI;
use glam::DVec3;
use lambert_izzo::{
    LambertSolution, RevolutionBudget, TransferWay, lambert, solve_with_diagnostics,
};
use lambert_izzo_test_support::bodies::{AU, MU_EARTH, MU_SUN};

fn print_trajectory(label: &str, s: LambertSolution) {
    let v1 = DVec3::from_array(s.v1);
    let v2 = DVec3::from_array(s.v2);
    println!(
        "  {label}: v1 = [{:+.6}, {:+.6}, {:+.6}] km/s  |v1| = {:.4} km/s  |v2| = {:.4} km/s",
        v1.x, v1.y, v1.z, v1.length(), v2.length(),
    );
}

fn main() {
    // 1. LEO → MEO Hohmann transfer (Earth-centered).
    println!("=== LEO (7000 km) → MEO (12000 km) Hohmann, short way ===");
    let r1 = [7000.0, 0.0, 0.0];
    let r2 = [-12_000.0, 1.0, 0.0]; // 1 km off-axis for non-collinearity
    let a = f64::midpoint(7000.0, 12_000.0);
    let tof = PI * (a.powi(3) / MU_EARTH).sqrt();
    let sols = lambert(
        r1,
        r2,
        tof,
        MU_EARTH,
        TransferWay::Short,
        RevolutionBudget::SingleOnly,
    )
    .expect("LEO Hohmann should converge");
    print_trajectory("single-rev", sols.single);

    // 2. Same geometry, long way around.
    println!("\n=== Same geometry, long way ===");
    let sols = lambert(
        r1,
        r2,
        tof,
        MU_EARTH,
        TransferWay::Long,
        RevolutionBudget::SingleOnly,
    )
    .expect("long-way Hohmann should converge");
    print_trajectory("single-rev", sols.single);

    // 3. Earth → Mars heliocentric Hohmann (Sun-centered, large scale).
    println!("\n=== Earth (1 AU) → Mars (1.524 AU) heliocentric, short way ===");
    let r1 = [AU, 0.0, 0.0];
    let r2 = [-1.524 * AU, 1.0, 0.0];
    let a = f64::midpoint(AU, 1.524 * AU);
    let tof = PI * (a.powi(3) / MU_SUN).sqrt();
    let sols = lambert(
        r1,
        r2,
        tof,
        MU_SUN,
        TransferWay::Short,
        RevolutionBudget::SingleOnly,
    )
    .expect("Earth-Mars Hohmann should converge");
    print_trajectory(
        &format!("single-rev (tof = {:.2} days)", tof / 86400.0),
        sols.single,
    );

    // 4. Earth-orbit multi-rev phasing — show the new pair structure.
    println!("\n=== Earth-orbit phasing (long tof, M up to 3) ===");
    let r1 = [8000.0, 0.0, 0.0];
    let r2 = [5600.0, 5600.0, 0.0];
    let period = 2.0 * PI * (8000.0_f64.powi(3) / MU_EARTH).sqrt();
    let tof = 5.0 * period;
    let (sols, diag) = solve_with_diagnostics(
        r1,
        r2,
        tof,
        MU_EARTH,
        TransferWay::Short,
        RevolutionBudget::up_to(3),
    )
    .expect("multi-rev phasing should converge");
    println!(
        "  single-rev: iters={} x={:+.6}  |v1|={:.4} km/s",
        diag.single.iters,
        diag.single.lancaster_blanchard_x,
        DVec3::from_array(sols.single.v1).length(),
    );
    for (pair, dpair) in sols.multi.iter().zip(diag.multi.iter()) {
        println!(
            "  M={} long-period:  iters={} x={:+.6}  |v1|={:.4} km/s",
            pair.n_revs,
            dpair.long_period.iters,
            dpair.long_period.lancaster_blanchard_x,
            DVec3::from_array(pair.long_period.v1).length(),
        );
        println!(
            "  M={} short-period: iters={} x={:+.6}  |v1|={:.4} km/s",
            pair.n_revs,
            dpair.short_period.iters,
            dpair.short_period.lancaster_blanchard_x,
            DVec3::from_array(pair.short_period.v1).length(),
        );
    }
}
