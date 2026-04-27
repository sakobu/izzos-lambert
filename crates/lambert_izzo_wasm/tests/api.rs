use lambert_izzo_wasm::{LambertRequest, TransferWayInput, solve_lambert_request};

#[test]
fn single_rev_request_returns_js_friendly_vectors() {
    let request = LambertRequest {
        r1_km: [7000.0, 0.0, 0.0],
        r2_km: [0.0, 7000.0, 0.0],
        tof_s: core::f64::consts::PI / 2.0 * (7000.0_f64.powi(3) / 398_600.441_8).sqrt(),
        mu_km3_s2: 398_600.441_8,
        way: TransferWayInput::Short,
        max_revs: 0,
    };

    let response = solve_lambert_request(request).unwrap();

    assert_eq!(response.solutions.len(), 1);
    let solution = &response.solutions[0];
    assert_eq!(solution.n_revs, 0);
    assert!(
        solution
            .v1_km_s
            .iter()
            .all(|component| component.is_finite())
    );
    assert!(
        solution
            .v2_km_s
            .iter()
            .all(|component| component.is_finite())
    );
    assert!(solution.diagnostics.iters > 0);
    assert!(solution.diagnostics.x.is_finite());
}

#[test]
fn invalid_request_returns_error_string() {
    let request = LambertRequest {
        r1_km: [0.0, 0.0, 0.0],
        r2_km: [0.0, 7000.0, 0.0],
        tof_s: 1000.0,
        mu_km3_s2: 398_600.441_8,
        way: TransferWayInput::Short,
        max_revs: 0,
    };

    let error = solve_lambert_request(request).unwrap_err();

    assert!(error.contains("degenerate position vector r1"));
}
