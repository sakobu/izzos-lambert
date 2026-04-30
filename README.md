# lambert_izzo

Cargo workspace for Izzo's revisited Lambert solver and its WebAssembly
adapter. A standalone Lambert crate — designed for callers who want a
focused, correct, `no_std`/WASM-friendly solver without pulling a full
astrodynamics framework just for the boundary-value step.

Reference: D. Izzo, *Revisiting Lambert's problem*, Celestial Mechanics &
Dynamical Astronomy, 2014 (arXiv:1403.2705).

MSRV: **Rust 1.85** (first edition-2024 stable). `rust-toolchain.toml`
pins **1.88.0** for development; the published-crate MSRV is 1.85.

## Crates

| Crate | Path | Distribution |
|-------|------|--------------|
| `lambert_izzo` | [`crates/lambert_izzo`](crates/lambert_izzo/README.md) | crates.io — pure Rust, `[f64; 3]` API, no hard math dep, `no_std`. |
| `lambert_izzo_wasm` | [`crates/lambert_izzo_wasm`](crates/lambert_izzo_wasm/README.md) | npm via `wasm-pack` — JS/TS bindings. Published as [`lambert-izzo`](https://www.npmjs.com/package/lambert-izzo). |
| `lambert_izzo_test_support` | [`crates/lambert_izzo_test_support`](crates/lambert_izzo_test_support/README.md) | workspace-internal (`publish = false`) — astrodynamics constants, Kepler propagator, deterministic random batches for examples / benches / tests. |

The core crate stays free of JavaScript concerns. The WASM crate owns
TypeScript type generation and JavaScript error conversion. The test
support crate is path-only — it never reaches `crates.io`.

See [`ARCHITECTURE.md`](ARCHITECTURE.md) for a code map of the solver
pipeline (input validation → geometry → root finding → TOF dispatch →
velocity reconstruction).

## Rust usage

```rust
use lambert_izzo::{lambert, LambertInput, RevolutionBudget, TransferWay};

let input = LambertInput {
    r1: [7000.0, 0.0, 0.0],
    r2: [0.0, 7000.0, 0.0],
    tof: 1457.0, // ~quarter-period of a 7000 km circular Earth orbit
    mu: 398_600.4418,
    way: TransferWay::Short,
    revolutions: RevolutionBudget::SingleOnly,
};

let solutions = lambert(&input)?;
let v1 = solutions.single.v1;
let iters = solutions.diagnostics.single.iters;
# Ok::<(), lambert_izzo::LambertError>(())
```

For batch / porkchop-plot workloads, enable the `rayon` feature and pass a
slice through `lambert_par`.

## WASM usage

Build the wrapper as an npm package with `wasm-pack`:

```bash
wasm-pack build crates/lambert_izzo_wasm --target bundler --release
```

`--target bundler` matches the published npm package and works with
Vite, webpack, Rollup, esbuild, and any other bundler that resolves
`import` statements. For direct browser ES-module imports without a
bundler, swap in `--target web`; for CommonJS Node.js consumption, use
`--target nodejs`. The full target table is in
[`crates/lambert_izzo_wasm/README.md`](crates/lambert_izzo_wasm/README.md).

Then import the generated package from a browser or bundler app:

```ts
import init, { solveLambert } from "./pkg/lambert_izzo_wasm";

await init();

const response = solveLambert({
  r1: [7000, 0, 0],
  r2: [0, 7000, 0],
  tof: 1457, // ~quarter-period of a 7000 km circular Earth orbit
  mu: 398600.4418,
  way: "short",
  maxRevs: null, // null = single-rev only; pass 1..=32 to search multi-rev branches
});

console.log(response.single.v1);
console.log(response.diagnostics.single.iters);
```

`maxRevs` is validated at request time. Out-of-range values (`0` or `> 32`)
reject with a typed `RevsOutOfRange` error before any solver work runs:

```ts
try {
  solveLambert({ ...request, maxRevs: 100 });
} catch (error) {
  // error.kind === "RevsOutOfRange"
  // error.requested === 100, error.max === 32
}
```

The wrapper response shape (camelCased per `serde(rename_all)`):

```ts
type LambertResponse = {
  single: { v1: [number, number, number]; v2: [number, number, number] };
  multi: Array<{
    nRevs: number;
    longPeriod: { v1: [...]; v2: [...] };
    shortPeriod: { v1: [...]; v2: [...] };
  }>;
  diagnostics: {
    single: { iters: number };
    multi: Array<{ nRevs: number; longPeriod: { iters: number }; shortPeriod: { iters: number } }>;
  };
};
```

A single-page browser demo lives at `crates/lambert_izzo_wasm/examples/web/`.

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

## Release notes

See [`CHANGELOG.md`](CHANGELOG.md) for per-release changes.

## License

MIT OR Apache-2.0
