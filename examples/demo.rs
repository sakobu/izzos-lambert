//! Worked examples for the [`lambert_izzo`] crate.
//!
//! Run with `cargo run --release --example demo`.

use core::f64::consts::PI;
use lambert_izzo::{lambert, RevolutionBudget, TransferWay};
use nalgebra::Vector3;

const MU_EARTH_KM3_S2: f64 = 398_600.441_8;
const MU_SUN_KM3_S2: f64 = 1.327_124_400_18e11;
const AU_KM: f64 = 1.495_978_707e8;

fn main() {
    // 1. LEO → MEO Hohmann transfer (Earth-centered).
    println!("=== LEO (7000 km) → MEO (12000 km) Hohmann, short way ===");
    let r1_km = Vector3::new(7000.0, 0.0, 0.0);
    let r2_km = Vector3::new(-12_000.0, 1.0, 0.0); // 1 km off-axis for non-collinearity
    let a_km = f64::midpoint(7000.0, 12_000.0);
    let tof_s = PI * (a_km.powi(3) / MU_EARTH_KM3_S2).sqrt();
    for s in lambert(r1_km, r2_km, tof_s, MU_EARTH_KM3_S2, TransferWay::Short, RevolutionBudget::SingleOnly).unwrap() {
        println!(
            "  M={} iters={} x={:+.6}\n    v1_km_s = [{:+.6}, {:+.6}, {:+.6}]\n    v2_km_s = [{:+.6}, {:+.6}, {:+.6}]\n    |v1|={:.4} km/s  |v2|={:.4} km/s",
            s.n_revs, s.iters, s.x,
            s.v1_km_s.x, s.v1_km_s.y, s.v1_km_s.z,
            s.v2_km_s.x, s.v2_km_s.y, s.v2_km_s.z,
            s.v1_km_s.norm(), s.v2_km_s.norm(),
        );
    }

    // 2. Same geometry, long way around.
    println!("\n=== Same geometry, long way ===");
    for s in lambert(r1_km, r2_km, tof_s, MU_EARTH_KM3_S2, TransferWay::Long, RevolutionBudget::SingleOnly).unwrap() {
        println!(
            "  M={} iters={} x={:+.6}  |v1|={:.4} km/s  |v2|={:.4} km/s",
            s.n_revs, s.iters, s.x, s.v1_km_s.norm(), s.v2_km_s.norm(),
        );
    }

    // 3. Earth → Mars heliocentric Hohmann (Sun-centered, large scale).
    println!("\n=== Earth (1 AU) → Mars (1.524 AU) heliocentric, short way ===");
    let r1_km = Vector3::new(AU_KM, 0.0, 0.0);
    let r2_km = Vector3::new(-1.524 * AU_KM, 1.0, 0.0);
    let a_km = f64::midpoint(AU_KM, 1.524 * AU_KM);
    let tof_s = PI * (a_km.powi(3) / MU_SUN_KM3_S2).sqrt();
    for s in lambert(r1_km, r2_km, tof_s, MU_SUN_KM3_S2, TransferWay::Short, RevolutionBudget::SingleOnly).unwrap() {
        println!(
            "  M={} iters={} x={:+.6}  |v1|={:.4} km/s  |v2|={:.4} km/s  (tof = {:.2} days)",
            s.n_revs, s.iters, s.x, s.v1_km_s.norm(), s.v2_km_s.norm(),
            tof_s / 86400.0,
        );
    }

    // 4. Earth-orbit multi-rev phasing.
    println!("\n=== Earth-orbit phasing (long tof, M up to 3) ===");
    let r1_km = Vector3::new(8000.0, 0.0, 0.0);
    let r2_km = Vector3::new(5600.0, 5600.0, 0.0);
    let period_s = 2.0 * PI * (8000.0_f64.powi(3) / MU_EARTH_KM3_S2).sqrt();
    let tof_s = 5.0 * period_s;
    for s in lambert(r1_km, r2_km, tof_s, MU_EARTH_KM3_S2, TransferWay::Short, RevolutionBudget::up_to(3)).unwrap() {
        println!(
            "  M={} iters={} x={:+.6}  |v1|={:.4} km/s  |v2|={:.4} km/s",
            s.n_revs, s.iters, s.x, s.v1_km_s.norm(), s.v2_km_s.norm(),
        );
    }
}
