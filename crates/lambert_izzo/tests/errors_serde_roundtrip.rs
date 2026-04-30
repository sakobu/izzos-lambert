//! Wire-format regression test for `LambertError`.
//!
//! Asserts every variant of the public `LambertError` enum round-trips
//! losslessly through `serde_json`. Catches accidental schema drift —
//! variant renames, field-name changes, or custom `#[serde(...)]`
//! attributes that would break existing JSON consumers. Lives at
//! `tests/` (integration-test boundary) so it also compile-checks the
//! public re-exports it imports.

#![cfg(feature = "serde")]

use lambert_izzo::{LambertError, NonFiniteParameter, Position};

#[test]
fn errors_round_trip_through_serde_json() {
    // `serde_json` rejects `f64::NAN` and `f64::INFINITY` outright
    // (standard JSON has no representation for them), so the
    // `NonFiniteInput` case uses `f64::MAX`. The variant's contract is
    // "non-finite at solver entry"; for serde round-trip purposes any
    // finite payload exercises the wire format identically.
    let cases = [
        LambertError::NonFiniteInput {
            parameter: NonFiniteParameter::Tof,
            value: f64::MAX,
        },
        LambertError::NonPositiveTimeOfFlight { tof: -1.0 },
        LambertError::NonPositiveMu { mu: 0.0 },
        LambertError::DegeneratePositionVector {
            position: Position::R1,
            norm: 0.0,
        },
        LambertError::CollinearGeometry { sin_angle: 1e-16 },
        LambertError::NoConvergence {
            iterations: 15,
            last_step: 1e-3,
            n_revs: 0,
        },
        LambertError::NoConvergence {
            iterations: 15,
            last_step: 1e-7,
            n_revs: 4,
        },
        LambertError::SingularDenominator { n_revs: 0 },
        LambertError::SingularDenominator { n_revs: 5 },
    ];

    for original in &cases {
        let json = serde_json::to_string(original)
            .expect("LambertError serializes to JSON");
        let restored: LambertError = serde_json::from_str(&json)
            .expect("LambertError deserializes from its own JSON");
        assert_eq!(original, &restored, "round-trip mismatch via {json}");
    }
}
