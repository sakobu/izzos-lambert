# lambert_izzo

Cargo workspace for Izzo's revisited Lambert solver and its WebAssembly
adapter.

## Crates

- `crates/lambert_izzo`: pure Rust solver crate. This is the canonical Rust
  API and the crate intended for normal Rust consumers. Public surface is
  `[f64; 3]` arrays — no hard math-library dependency.
- `crates/lambert_izzo_wasm`: thin `wasm-bindgen` adapter that exposes
  JavaScript and TypeScript friendly request and response types.

The core crate stays free of JavaScript concerns. The WASM crate owns
TypeScript type generation and JavaScript error conversion.

## Rust usage

```rust
use lambert_izzo::{lambert, RevolutionBudget, TransferWay};

let mu_km3_s2 = 398_600.441_8;
let r1_km = [7000.0, 0.0, 0.0];
let r2_km = [0.0, 7000.0, 0.0];
let tof_s = core::f64::consts::PI / 2.0 * (7000.0_f64.powi(3) / mu_km3_s2).sqrt();

let solutions = lambert(
    r1_km,
    r2_km,
    tof_s,
    mu_km3_s2,
    TransferWay::Short,
    RevolutionBudget::SingleOnly,
)?;
let v1_km_s = solutions.single.v1_km_s;
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
  r1Km: [7000, 0, 0],
  r2Km: [0, 7000, 0],
  tofS: 1457,
  muKm3S2: 398600.4418,
  way: "short",
  maxRevs: 0,
});

console.log(response.single.v1KmS);
```

The wrapper returns the same shape as the Rust core (camelCased):

```ts
{
  single: {
    v1KmS: [number, number, number],
    v2KmS: [number, number, number],
    diagnostics: { iters: number, x: number }
  },
  multi: [
    {
      nRevs: number,
      longPeriod:  { v1KmS, v2KmS, diagnostics },
      shortPeriod: { v1KmS, v2KmS, diagnostics }
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
