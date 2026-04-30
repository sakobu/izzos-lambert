//! Multi-revolution Lambert throughput (`M = 3`).
//!
//! Run: `cargo bench --bench multi_rev`.

use criterion::{Criterion, criterion_group, criterion_main};
use lambert_izzo::{RevolutionBudget, lambert};
use lambert_izzo_test_support::bodies::MU_EARTH;
use lambert_izzo_test_support::random_inputs::{Spec, WayStrategy, generate};
use std::hint::black_box;

fn multi_rev_throughput(c: &mut Criterion) {
    let inputs = generate(&Spec {
        n: 1_000,
        seed: 0xBEEF_DEAD,
        radius_range: (5600.0, 10_500.0),
        tof_range: (10_000.0, 250_000.0),
        mu: MU_EARTH,
        way: WayStrategy::Short,
        revolutions: RevolutionBudget::up_to(3),
    });

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
