# lambert_izzo

Cargo workspace for Izzo's revisited Lambert solver and its WebAssembly
adapter. A standalone Lambert crate — designed for callers who want a
focused, correct, `no_std`/WASM-friendly solver without pulling a full
astrodynamics framework just for the boundary-value step.

MSRV: **Rust 1.85** (first edition-2024 stable). Pre-1.0; breaking
changes are tracked in [`CHANGELOG.md`](CHANGELOG.md).

## Crates

- `crates/lambert_izzo`: pure Rust solver crate. This is the canonical Rust
  API and the crate intended for normal Rust consumers. Public surface is
  `[f64; 3]` arrays — no hard math-library dependency.
- `crates/lambert_izzo_wasm`: thin `wasm-bindgen` adapter that exposes
  JavaScript and TypeScript friendly request and response types.
- `crates/lambert_izzo_test_support`: workspace-internal dev fixtures
  (`publish = false`) — astrodynamics constants (`MU_EARTH`, `MU_SUN`,
  `AU`), a rejection-sampling unit-vector helper, and a
  universal-variable Kepler propagator used by the examples, benches,
  and integration tests across both other crates.

The core crate stays free of JavaScript concerns. The WASM crate owns
TypeScript type generation and JavaScript error conversion. The test
support crate is path-only — it never reaches `crates.io`.

## Rust usage

```rust
use lambert_izzo::{lambert, RevolutionBudget, TransferWay};

let mu = 398_600.441_8;
let r1 = [7000.0, 0.0, 0.0];
let r2 = [0.0, 7000.0, 0.0];
let tof = core::f64::consts::PI / 2.0 * (7000.0_f64.powi(3) / mu).sqrt();

let solutions = lambert(
    r1,
    r2,
    tof,
    mu,
    TransferWay::Short,
    RevolutionBudget::SingleOnly,
)?;
let v1 = solutions.single.v1;
```

## WASM usage

Build the wrapper as an npm package with `wasm-pack`:

```bash
wasm-pack build crates/lambert_izzo_wasm --target bundler --out-dir ../../pkg
```

Then import the generated package from a browser or bundler app:

```ts
import init, { solveLambert } from "./pkg/lambert_izzo_wasm";

await init();

const response = solveLambert({
  r1: [7000, 0, 0],
  r2: [0, 7000, 0],
  tof: 1457,
  mu: 398600.4418,
  way: "short",
  maxRevs: 0,
});

console.log(response.single.v1);
```

The wrapper returns the same shape as the Rust core (camelCased):

```ts
{
  single: {
    v1: [number, number, number],
    v2: [number, number, number],
    diagnostics: { iters: number, x: number }
  },
  multi: [
    {
      nRevs: number,
      longPeriod:  { v1, v2, diagnostics },
      shortPeriod: { v1, v2, diagnostics }
    }
  ]
}
```

## Development

```bash
cargo build --workspace
cargo test --workspace
cargo clippy --workspace --all-targets -- -D warnings
cargo build --target wasm32-unknown-unknown -p lambert_izzo --lib
cargo build --target wasm32-unknown-unknown -p lambert_izzo_wasm --lib
cargo run --release -p lambert_izzo --example demo
cargo run --release -p lambert_izzo --example stress
```

Toolchain is pinned by `rust-toolchain.toml`.

## License

MIT OR Apache-2.0
