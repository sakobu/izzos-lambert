//! Single-revolution Lambert throughput.
//!
//! Run: `cargo bench --bench single_rev`.

use criterion::{Criterion, criterion_group, criterion_main};
use lambert_izzo::{RevolutionBudget, lambert};
use lambert_izzo_test_support::bodies::MU_EARTH;
use lambert_izzo_test_support::random_inputs::{Spec, WayStrategy, generate};
use std::hint::black_box;

fn single_rev_throughput(c: &mut Criterion) {
    let inputs = generate(&Spec {
        n: 10_000,
        seed: 0xC0FF_EE42,
        radius_range: (3500.0, 28_000.0),
        tof_range: (100.0, 50_000.0),
        mu: MU_EARTH,
        way: WayStrategy::Random,
        revolutions: RevolutionBudget::SingleOnly,
    });

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
