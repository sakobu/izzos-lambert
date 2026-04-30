# lambert_izzo_test_support

Workspace-internal dev fixtures for `lambert_izzo`. **Not published**
(`publish = false`) and **not stable** — no semver guarantees on its
public surface. Used as a `[dev-dependencies]` path crate by the core's
examples, benches, and tests.

## What's in it

- **`bodies`** — SI astrodynamics constants: `MU_EARTH`, `MU_SUN`, `AU`.
- **`rand_unit_vec`** — uniform-on-sphere unit-vector helper via
  rejection sampling on the unit cube.
- **`random_inputs`** — deterministic `LambertInput` batch generator
  (`Spec`, `generate`, `WayStrategy`). Seeded with `ChaCha20Rng` so
  benches and stress runs are reproducible.
- **`kepler::propagate`** — universal-variable / Stumpff Kepler
  propagator used to round-trip Lambert solutions in the stress example
  and in `tests/kepler_roundtrip.rs`.
- **`elements_u64`** — `usize → u64` converter for
  `criterion::Throughput::Elements`, casting policy-compliant.

## Why a separate crate?

Keeps `lambert_izzo` itself dependency-free at runtime: nothing here
(`rand`, `rand_chacha`, the SI constants) leaks into the published
crate. The core can stay `no_std` while the examples / benches use the
real `std` test scaffolding.

## Bar for changes

Anything added here should pass the bar for `lambert_izzo` proper —
casting policy, no `unwrap`/`expect`/`panic!` outside test assertions,
explanatory `///` doc comments on every public item — even though this
crate isn't published.
