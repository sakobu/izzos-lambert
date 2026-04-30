//! Worked tour of `LambertError` — when each variant fires and how to
//! recover. Run with `cargo run --release --example errors`.

use core::f64::consts::PI;

use glam::DVec3;
use lambert_izzo::{
    LambertError, LambertInput, NonFiniteParameter, Position, RevolutionBudget, TransferWay,
    lambert,
};
use lambert_izzo_test_support::bodies::MU_EARTH;

fn main() {
    println!("=== 1. Colinear geometry — perturb to recover ===");
    // r1 and r2 both on the +X axis: transfer plane undefined.
    let r1 = [7000.0, 0.0, 0.0];
    let r2 = [14_000.0, 0.0, 0.0];
    match lambert(&LambertInput {
        r1,
        r2,
        tof: 1500.0,
        mu: MU_EARTH,
        way: TransferWay::Short,
        revolutions: RevolutionBudget::SingleOnly,
    }) {
        Err(LambertError::CollinearGeometry { sin_angle }) => {
            println!(
                "  rejected: sin(θ) = {sin_angle:.3e} (transfer plane undefined).\n  \
                 Recovery: perturb r2 by 1 km off-plane."
            );
            // Re-solve with a tiny off-plane component.
            let r2_perturbed = [r2[0], 1.0, 0.0];
            let sols = lambert(&LambertInput {
                r1,
                r2: r2_perturbed,
                tof: 1500.0,
                mu: MU_EARTH,
                way: TransferWay::Short,
                revolutions: RevolutionBudget::SingleOnly,
            })
            .expect("perturbed geometry should converge");
            println!(
                "  perturbed v1 = [{:+.4}, {:+.4}, {:+.4}] km/s",
                sols.single.v1[0], sols.single.v1[1], sols.single.v1[2],
            );
        }
        Err(e) => println!("  unexpected error: {e}"),
        Ok(_) => println!("  unexpected success"),
    }

    println!("\n=== 2. Non-finite input — typed Position + parameter ===");
    let bad_r1 = [7000.0, f64::INFINITY, 0.0];
    let r2 = [0.0, 7000.0, 0.0];
    match lambert(&LambertInput {
        r1: bad_r1,
        r2,
        tof: 1500.0,
        mu: MU_EARTH,
        way: TransferWay::Short,
        revolutions: RevolutionBudget::SingleOnly,
    }) {
        Err(LambertError::NonFiniteInput { parameter, value }) => {
            // Pattern-match on the typed enum, not a string.
            let component = match parameter {
                NonFiniteParameter::R1Y => "r1.y",
                other => other.as_str(),
            };
            println!("  rejected: {component} = {value} (must be finite).");
        }
        other => println!("  unexpected: {other:?}"),
    }

    println!("\n=== 3. Multi-rev infeasibility — silent skip ===");
    // Earth-orbit phasing with up_to(10), but tof only large enough for a
    // few branches. Higher M get dropped from the returned `multi` set.
    let r1 = [8000.0, 0.0, 0.0];
    let r2 = [5600.0, 5600.0, 0.0];
    let period = 2.0 * PI * (8000.0_f64.powi(3) / MU_EARTH).sqrt();
    let tof = 3.0 * period;
    let sols = lambert(&LambertInput {
        r1,
        r2,
        tof,
        mu: MU_EARTH,
        way: TransferWay::Short,
        revolutions: RevolutionBudget::up_to(10),
    })
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
            DVec3::from_array(pair.long_period.v1).length(),
            DVec3::from_array(pair.short_period.v1).length(),
        );
    }

    println!("\n=== 4. Near-parabolic — Battin dispatch is automatic ===");
    let r1 = [7000.0, 0.0, 0.0];
    let r2 = [0.0, 42_000.0, 0.0];
    let tof = 7200.0;
    let sols = lambert(&LambertInput {
        r1,
        r2,
        tof,
        mu: MU_EARTH,
        way: TransferWay::Short,
        revolutions: RevolutionBudget::SingleOnly,
    })
    .expect("near-parabolic should converge via Battin");
    println!(
        "  iters = {} (single-rev paper avg ≈ 2.1); regime is selected internally \
         when |x − 1| is small enough to need the hypergeometric series.",
        sols.diagnostics.single.iters,
    );
    println!(
        "  trajectory v1 = [{:+.4}, {:+.4}, {:+.4}] km/s",
        sols.single.v1[0], sols.single.v1[1], sols.single.v1[2],
    );

    println!("\n=== 5. Degenerate position — typed Position::R1 / R2 ===");
    let err = lambert(&LambertInput {
        r1: [0.0, 0.0, 0.0],
        r2,
        tof: 1500.0,
        mu: MU_EARTH,
        way: TransferWay::Short,
        revolutions: RevolutionBudget::SingleOnly,
    })
    .unwrap_err();
    if let LambertError::DegeneratePositionVector { position, norm } = err {
        let which = match position {
            Position::R1 => "r1",
            Position::R2 => "r2",
        };
        println!("  rejected: {which} has norm {norm:.3e} (below internal floor).");
    }

    println!("\n=== 6. Non-positive scalar — distinct error variants ===");
    for (label, err) in [
        (
            "tof = 0",
            lambert(&LambertInput {
                r1,
                r2,
                tof: 0.0,
                mu: MU_EARTH,
                way: TransferWay::Short,
                revolutions: RevolutionBudget::SingleOnly,
            }),
        ),
        (
            "mu = -1",
            lambert(&LambertInput {
                r1,
                r2,
                tof: 1500.0,
                mu: -1.0,
                way: TransferWay::Short,
                revolutions: RevolutionBudget::SingleOnly,
            }),
        ),
    ] {
        match err {
            Err(LambertError::NonPositiveTimeOfFlight { tof }) => {
                println!("  {label}: NonPositiveTimeOfFlight(tof = {tof})");
            }
            Err(LambertError::NonPositiveMu { mu }) => {
                println!("  {label}: NonPositiveMu(mu = {mu})");
            }
            Err(other) => println!("  {label}: unexpected variant {other:?}"),
            Ok(_) => println!("  {label}: unexpected success"),
        }
    }
}
