//! Batch throughput: sequential `lambert_iter` vs (with `--features rayon`)
//! parallel `lambert_par_iter`.
//!
//! Run sequential: `cargo bench --bench batch`.
//! Run parallel:   `cargo bench --bench batch --features rayon`.

use criterion::{Criterion, criterion_group, criterion_main};
use lambert_izzo::{LambertInput, RevolutionBudget, TransferWay, lambert_iter};
use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha20Rng;
use rand_distr::Uniform;
use std::hint::black_box;

const MU_EARTH: f64 = 398_600.441_8;

fn rand_unit_vec(rng: &mut ChaCha20Rng) -> [f64; 3] {
    let axis: Uniform<f64> = Uniform::new(-1.0, 1.0);
    loop {
        let v: [f64; 3] = [rng.sample(axis), rng.sample(axis), rng.sample(axis)];
        let n2 = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
        if n2 > 0.01 && n2 < 1.0 {
            let inv = 1.0 / n2.sqrt();
            return [v[0] * inv, v[1] * inv, v[2] * inv];
        }
    }
}

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
                way: TransferWay::Short,
                revolutions: RevolutionBudget::up_to(3),
            }
        })
        .collect()
}

fn batch_sequential(c: &mut Criterion) {
    let inputs = build_inputs(10_000, 0xCAFE_F00D);

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

    let inputs = build_inputs(10_000, 0xCAFE_F00D);

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
