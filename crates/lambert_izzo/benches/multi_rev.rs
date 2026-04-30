//! Multi-revolution Lambert throughput across `M ∈ {1, 3, 5}` plus a
//! Battin near-parabolic regime stress.
//!
//! Two benchmark groups:
//! - `multi_rev`: random Earth-orbit inputs parametrized over the
//!   revolution budget via Criterion's `bench_with_input`.
//! - `multi_rev_battin`: a deterministic near-180° geometry whose
//!   solutions sit in the Battin regime (`|x − 1| < BATTIN_THRESHOLD`),
//!   constructed by jittering TOF around the parabolic time-of-flight.
//!
//! Run: `cargo bench --bench multi_rev`.

use criterion::{BenchmarkId, Criterion, criterion_group, criterion_main};
use lambert_izzo::{LambertInput, RevolutionBudget, TransferWay, lambert};
use lambert_izzo_test_support::bodies::MU_EARTH;
use lambert_izzo_test_support::elements_u64;
use lambert_izzo_test_support::random_inputs::{Spec, WayStrategy, generate};
use rand::SeedableRng;
use rand::distributions::{Distribution, Uniform};
use rand_chacha::ChaCha20Rng;
use std::hint::black_box;

fn multi_rev_throughput(c: &mut Criterion) {
    let mut group = c.benchmark_group("multi_rev");
    group.sample_size(20);

    for &m in &[1u32, 3, 5] {
        let inputs = generate(&Spec {
            n: 1_000,
            seed: 0xBEEF_DEAD ^ u64::from(m),
            radius_range: (5_600.0, 10_500.0),
            tof_range: (10_000.0, 250_000.0),
            mu: MU_EARTH,
            way: WayStrategy::Short,
            revolutions: RevolutionBudget::try_up_to(m).expect("m ≤ BoundedRevs::MAX"),
        });

        group.throughput(criterion::Throughput::Elements(elements_u64(inputs.len())));
        group.bench_with_input(BenchmarkId::new("M", m), &inputs, |b, inputs| {
            b.iter(|| {
                for input in inputs {
                    let _ = black_box(lambert(input));
                }
            });
        });
    }
    group.finish();
}

fn battin_near_parabolic(c: &mut Criterion) {
    let inputs = generate_battin_inputs(1_000, MU_EARTH);

    let mut group = c.benchmark_group("multi_rev_battin");
    group.throughput(criterion::Throughput::Elements(elements_u64(inputs.len())));
    group.sample_size(20);
    group.bench_function("near_parabolic_x1000", |b| {
        b.iter(|| {
            for input in &inputs {
                let _ = black_box(lambert(input));
            }
        });
    });
    group.finish();
}

/// Deterministic batch whose solutions land in `tof.rs`'s Battin regime
/// (`|x − 1| < BATTIN_THRESHOLD = 0.01`). Achieved via a fixed near-180°
/// transfer geometry with TOF jittered around the parabolic time-of-flight
/// (Lancaster's closed form `T_p = (√2 / (3√μ)) · (s^(3/2) − (s−c)^(3/2))`
/// for the short-way branch).
fn generate_battin_inputs(n: usize, mu: f64) -> Vec<LambertInput> {
    let r_norm = 7_000.0_f64;
    // ~177° transfer — close enough to parabolic to land in Battin,
    // far enough from 180° to clear `COLINEARITY_TOL = 1e-15`.
    let theta = std::f64::consts::PI - 0.05;
    let r1 = [r_norm, 0.0, 0.0];
    let r2 = [r_norm * theta.cos(), r_norm * theta.sin(), 0.0];

    let dx = r2[0] - r1[0];
    let dy = r2[1] - r1[1];
    let chord = (dx * dx + dy * dy).sqrt();
    let s = (2.0 * r_norm + chord) / 2.0;
    let t_par = (2.0 / mu).sqrt() / 3.0 * (s.powf(1.5) - (s - chord).powf(1.5));

    let mut rng = ChaCha20Rng::seed_from_u64(0xBA77_171F);
    let dist = Uniform::new(0.995 * t_par, 1.005 * t_par);
    (0..n)
        .map(|_| LambertInput {
            r1,
            r2,
            tof: dist.sample(&mut rng),
            mu,
            way: TransferWay::Short,
            revolutions: RevolutionBudget::SingleOnly,
        })
        .collect()
}

criterion_group!(benches, multi_rev_throughput, battin_near_parabolic);
criterion_main!(benches);
