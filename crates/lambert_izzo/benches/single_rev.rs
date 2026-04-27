//! Single-revolution Lambert throughput.
//!
//! Run: `cargo bench --bench single_rev`.

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
    let radius = Uniform::new(3500.0_f64, 28_000.0_f64);
    let tof_dist = Uniform::new(100.0_f64, 50_000.0_f64);
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
                way: if rng.gen_bool(0.5) {
                    TransferWay::Long
                } else {
                    TransferWay::Short
                },
                revolutions: RevolutionBudget::SingleOnly,
            }
        })
        .collect()
}

fn single_rev_throughput(c: &mut Criterion) {
    let inputs = build_inputs(10_000, 0xC0FF_EE42);

    let mut group = c.benchmark_group("single_rev");
    group.throughput(criterion::Throughput::Elements(inputs.len() as u64));
    group.sample_size(20);
    group.bench_function("lambert_x10000_random_earth", |b| {
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

criterion_group!(benches, single_rev_throughput);
criterion_main!(benches);
