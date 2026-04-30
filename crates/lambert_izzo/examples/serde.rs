//! Serde round-trip example. Run with:
//! `cargo run --release --example serde --features serde`
//!
//! Demonstrates the JSON shapes of `LambertSolutions` and `LambertError`
//! under the `serde` feature, and round-trips both through `serde_json`
//! to verify lossless serialization.

use lambert_izzo::{
    LambertError, LambertInput, LambertSolutions, RevolutionBudget, TransferWay, lambert,
};
use lambert_izzo_test_support::bodies::MU_EARTH;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("=== 1. LambertSolutions happy-path round-trip ===");
    let input = LambertInput {
        r1: [7000.0, 0.0, 0.0],
        r2: [0.0, 20_000.0, 0.0],
        tof: 9500.0,
        mu: MU_EARTH,
        way: TransferWay::Short,
        revolutions: RevolutionBudget::SingleOnly,
    };
    let sols = lambert(&input).expect("LEO transfer should converge");

    let json = serde_json::to_string_pretty(&sols)?;
    println!("{json}");

    let restored: LambertSolutions = serde_json::from_str(&json)?;
    assert_eq!(sols, restored, "solutions must round-trip losslessly");
    println!("  round-trip OK (PartialEq)");

    println!("\n=== 2. LambertError discriminated-union round-trip ===");
    // NonPositiveTimeOfFlight carries a guaranteed-finite `tof: f64`; JSON
    // cannot represent NaN, so picking a NaN-payloaded variant here would
    // defeat the round-trip assertion below.
    let bad = LambertInput { tof: -1.0, ..input };
    let err = lambert(&bad).expect_err("negative tof must error");

    let err_json = serde_json::to_string_pretty(&err)?;
    println!("{err_json}");

    let restored_err: LambertError = serde_json::from_str(&err_json)?;
    assert_eq!(err, restored_err, "errors must round-trip losslessly");
    println!("  error round-trip OK (PartialEq)");

    Ok(())
}
