//! Batch throughput: sequential `lambert_iter` vs (with `--features rayon`)
//! parallel `lambert_par_iter`.
//!
//! Run sequential: `cargo bench --bench batch`.
//! Run parallel:   `cargo bench --bench batch --features rayon`.

use criterion::{Criterion, criterion_group, criterion_main};
use lambert_izzo::{RevolutionBudget, lambert_iter};
use lambert_izzo_test_support::bodies::MU_EARTH;
use lambert_izzo_test_support::random_inputs::{Spec, WayStrategy, generate};
use std::hint::black_box;

fn spec() -> Spec {
    Spec {
        n: 10_000,
        seed: 0xCAFE_F00D,
        radius_range: (3500.0, 28_000.0),
        tof_range: (100.0, 50_000.0),
        mu: MU_EARTH,
        way: WayStrategy::Short,
        revolutions: RevolutionBudget::up_to(3),
    }
}

fn batch_sequential(c: &mut Criterion) {
    let inputs = generate(&spec());

    let mut group = c.benchmark_group("batch_sequential");
    group.throughput(criterion::Throughput::Elements(inputs.len() as u64));
    group.sample_size(20);
    group.bench_function("lambert_iter_x10000", |b| {
        b.iter(|| {
            for sol in lambert_iter(&inputs) {
                black_box(sol.ok());
            }
        });
    });
    group.finish();
}

#[cfg(feature = "rayon")]
fn batch_parallel(c: &mut Criterion) {
    use lambert_izzo::lambert_par_iter;
    use rayon::iter::ParallelIterator;

    let inputs = generate(&spec());

    let mut group = c.benchmark_group("batch_parallel");
    group.throughput(criterion::Throughput::Elements(inputs.len() as u64));
    group.sample_size(20);
    group.bench_function("lambert_par_iter_x10000", |b| {
        b.iter(|| {
            lambert_par_iter(&inputs).for_each(|sol| {
                black_box(sol.ok());
            });
        });
    });
    group.finish();
}

#[cfg(feature = "rayon")]
criterion_group!(benches, batch_sequential, batch_parallel);
#[cfg(not(feature = "rayon"))]
criterion_group!(benches, batch_sequential);

criterion_main!(benches);
