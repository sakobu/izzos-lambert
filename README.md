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
