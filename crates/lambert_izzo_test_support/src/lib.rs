//! Workspace-internal dev fixtures used by `lambert_izzo`'s examples,
//! benches, and integration tests. Not published.

#![warn(missing_docs)]

use rand::Rng;

/// Standard astrodynamics constants in SI (km, km/s, s, km³/s²).
pub mod bodies {
    /// Earth's gravitational parameter (km³/s²) — value from EGM2008.
    pub const MU_EARTH: f64 = 398_600.441_8;

    /// Sun's gravitational parameter (km³/s²) — value from DE440.
    pub const MU_SUN: f64 = 1.327_124_400_18e11;

    /// One astronomical unit (km) — IAU 2012 definition.
    pub const AU: f64 = 1.495_978_707e8;
}

/// Generate a uniformly distributed unit vector by rejection sampling on
/// the unit cube. Used to seed Lambert stress / bench inputs.
pub fn rand_unit_vec<R: Rng + ?Sized>(rng: &mut R) -> [f64; 3] {
    loop {
        let x: f64 = rng.gen_range(-1.0..=1.0);
        let y: f64 = rng.gen_range(-1.0..=1.0);
        let z: f64 = rng.gen_range(-1.0..=1.0);
        let n2 = x * x + y * y + z * z;
        if n2 > 0.0 && n2 <= 1.0 {
            let n = n2.sqrt();
            return [x / n, y / n, z / n];
        }
    }
}

/// Deterministic random batches of [`lambert_izzo::LambertInput`] for
/// benches and stress runs.
///
/// Position vectors are uniformly distributed on the sphere (via
/// [`rand_unit_vec`]) and scaled by an independent uniform magnitude per
/// endpoint; `tof` is uniform on its range; `way` is configurable.
pub mod random_inputs {
    use lambert_izzo::{LambertInput, RevolutionBudget, TransferWay};
    use rand::distributions::Uniform;
    use rand::{Rng, SeedableRng};
    use rand_chacha::ChaCha20Rng;

    /// How to choose the transfer direction for each generated input.
    #[derive(Debug, Clone, Copy)]
    pub enum WayStrategy {
        /// Force [`TransferWay::Short`] for every input.
        Short,
        /// Force [`TransferWay::Long`] for every input.
        Long,
        /// 50/50 random Short/Long, drawn independently per input.
        Random,
    }

    /// Specification for a random `LambertInput` batch.
    ///
    /// Range tuples are `(min, max)` and follow `rand`'s `Uniform::new`
    /// semantics — inclusive of the lower bound, exclusive of the upper.
    #[derive(Debug, Clone, Copy)]
    pub struct Spec {
        /// How many inputs to produce.
        pub n: usize,
        /// RNG seed (deterministic batches across runs).
        pub seed: u64,
        /// `(min, max)` of the position-vector magnitude (km in SI).
        pub radius_range: (f64, f64),
        /// `(min, max)` of the time of flight (s in SI).
        pub tof_range: (f64, f64),
        /// Gravitational parameter (km³/s² in SI).
        pub mu: f64,
        /// How to pick the transfer way per input.
        pub way: WayStrategy,
        /// Revolution budget passed through to every input.
        pub revolutions: RevolutionBudget,
    }

    /// Build a deterministic batch of random `LambertInput`s following `spec`.
    #[must_use]
    pub fn generate(spec: &Spec) -> Vec<LambertInput> {
        let mut rng = ChaCha20Rng::seed_from_u64(spec.seed);
        let radius = Uniform::new(spec.radius_range.0, spec.radius_range.1);
        let tof = Uniform::new(spec.tof_range.0, spec.tof_range.1);
        (0..spec.n)
            .map(|_| {
                let r1n = rng.sample(radius);
                let r2n = rng.sample(radius);
                let r1u = crate::rand_unit_vec(&mut rng);
                let r2u = crate::rand_unit_vec(&mut rng);
                let way = match spec.way {
                    WayStrategy::Short => TransferWay::Short,
                    WayStrategy::Long => TransferWay::Long,
                    WayStrategy::Random => {
                        if rng.gen_bool(0.5) {
                            TransferWay::Long
                        } else {
                            TransferWay::Short
                        }
                    }
                };
                LambertInput {
                    r1: [r1u[0] * r1n, r1u[1] * r1n, r1u[2] * r1n],
                    r2: [r2u[0] * r2n, r2u[1] * r2n, r2u[2] * r2n],
                    tof: rng.sample(tof),
                    mu: spec.mu,
                    way,
                    revolutions: spec.revolutions,
                }
            })
            .collect()
    }
}

/// Reference Kepler propagator for Lambert round-trip validation.
pub mod kepler {
    /// Newton-iteration cap on the universal Kepler equation. Well above the
    /// `<30` typical convergence count; pathological geometries hit it and bail
    /// to NaN, which statistical tests filter explicitly.
    const MAX_ITERS: u32 = 200;

    /// Convergence tolerance on the universal-variable step `|Δχ|`.
    const NEWTON_TOL: f64 = 1e-12;

    /// Threshold below which the Newton denominator `r` is treated as zero.
    const DENOM_EPS: f64 = 1e-12;

    /// Divergence sentinel on `|χ|`. Caps universal-variable iteration so
    /// poor initial guesses bail to NaN rather than looping forever.
    const CHI_DIVERGENCE: f64 = 1e12;

    /// Branch threshold for the Stumpff `c2(ψ)`, `c3(ψ)` evaluator. Above
    /// this magnitude, use the closed-form trig/hyperbolic expressions;
    /// below it, use the Taylor series (closed forms suffer catastrophic
    /// cancellation as `ψ → 0`).
    const STUMPFF_SERIES_THRESHOLD: f64 = 1e-6;

    /// Propagate `(r0, v0)` for `dt` seconds under two-body gravity.
    ///
    /// Inputs in km, km/s, s, km³/s². Returns the propagated position in km,
    /// or a NaN-vector if the universal-variable Newton iteration diverges.
    /// Callers in statistical tests filter divergent trials by checking
    /// componentwise finiteness.
    #[must_use]
    pub fn propagate(r0: [f64; 3], v0: [f64; 3], dt: f64, mu: f64) -> [f64; 3] {
        let r0n = norm(r0);
        let v0n2 = dot(v0, v0);
        let alpha = 2.0 / r0n - v0n2 / mu;
        let sqrt_mu = mu.sqrt();
        let sigma0 = dot(r0, v0) / sqrt_mu;
        let mut chi = sqrt_mu * dt * alpha.abs();
        let nan_vec = [f64::NAN; 3];
        let mut converged = false;
        for _ in 0..MAX_ITERS {
            let psi = alpha * chi * chi;
            let (c2, c3) = stumpff(psi);
            let r =
                chi * chi * c2 + sigma0 * chi * (1.0 - psi * c3) + r0n * (1.0 - psi * c2);
            if r.abs() < DENOM_EPS || !r.is_finite() {
                return nan_vec;
            }
            let f = sigma0 * chi * chi * c2
                + (1.0 - alpha * r0n) * chi * chi * chi * c3
                + r0n * chi
                - sqrt_mu * dt;
            let dchi = -f / r;
            chi += dchi;
            if !chi.is_finite() || chi.abs() > CHI_DIVERGENCE {
                return nan_vec;
            }
            if dchi.abs() < NEWTON_TOL {
                converged = true;
                break;
            }
        }
        if !converged {
            return nan_vec;
        }
        let psi = alpha * chi * chi;
        let (c2, c3) = stumpff(psi);
        let f = 1.0 - chi * chi / r0n * c2;
        let g = dt - chi * chi * chi / sqrt_mu * c3;
        [
            r0[0] * f + v0[0] * g,
            r0[1] * f + v0[1] * g,
            r0[2] * f + v0[2] * g,
        ]
    }

    fn stumpff(psi: f64) -> (f64, f64) {
        if psi > STUMPFF_SERIES_THRESHOLD {
            let s = psi.sqrt();
            ((1.0 - s.cos()) / psi, (s - s.sin()) / s.powi(3))
        } else if psi < -STUMPFF_SERIES_THRESHOLD {
            let s = (-psi).sqrt();
            ((1.0 - s.cosh()) / psi, (s.sinh() - s) / s.powi(3))
        } else {
            // Taylor coefficients for c2(ψ) = (1 − cos√ψ)/ψ and c3(ψ) = (√ψ − sin√ψ)/√ψ³
            // expanded about ψ = 0; denominators are factorials (2!, 4!, 6!) and (3!, 5!, 7!).
            (
                0.5 - psi / 24.0 + psi * psi / 720.0,
                1.0 / 6.0 - psi / 120.0 + psi * psi / 5040.0,
            )
        }
    }

    fn norm(v: [f64; 3]) -> f64 {
        dot(v, v).sqrt()
    }

    fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
        a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
    }
}
