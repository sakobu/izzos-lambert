//! Multi-revolution Lambert throughput (`M = 3`).
//!
//! Run: `cargo bench --bench multi_rev`.

use criterion::{Criterion, criterion_group, criterion_main};
use lambert_izzo::{LambertInput, RevolutionBudget, TransferWay, lambert};
use lambert_izzo_test_support::bodies::MU_EARTH;
use lambert_izzo_test_support::rand_unit_vec;
use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha20Rng;
use rand_distr::Uniform;
use std::hint::black_box;

fn build_inputs(n: usize, seed: u64) -> Vec<LambertInput> {
    let mut rng = ChaCha20Rng::seed_from_u64(seed);
    let radius = Uniform::new(5600.0_f64, 10_500.0_f64);
    let tof_dist = Uniform::new(10_000.0_f64, 250_000.0_f64);
    (0..n)
        .map(|_| {
            let r1n = rng.sample(radius);
            let r2n = rng.sample(radius);
            let r1u = rand_unit_vec(&mut rng);
            let r2u = rand_unit_vec(&mut rng);
            LambertInput {
                r1: [r1u[0] * r1n, r1u[1] * r1n, r1u[2] * r1n],
                r2: [r2u[0] * r2n, r2u[1] * r2n, r2u[2] * r2n],
                tof: rng.sample(tof_dist),
                mu: MU_EARTH,
                way: TransferWay::Short,
                revolutions: RevolutionBudget::up_to(3),
            }
        })
        .collect()
}

fn multi_rev_throughput(c: &mut Criterion) {
    let inputs = build_inputs(1_000, 0xBEEF_DEAD);

    let mut group = c.benchmark_group("multi_rev_M3");
    group.throughput(criterion::Throughput::Elements(inputs.len() as u64));
    group.sample_size(20);
    group.bench_function("lambert_x1000_random_earth", |b| {
        b.iter(|| {
            for input in &inputs {
                let _ = black_box(lambert(
                    input.r1,
                    input.r2,
                    input.tof,
                    input.mu,
                    input.way,
                    input.revolutions,
                ));
            }
        });
    });
    group.finish();
}

criterion_group!(benches, multi_rev_throughput);
criterion_main!(benches);
