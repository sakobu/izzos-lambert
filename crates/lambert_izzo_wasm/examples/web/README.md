# lambert-izzo web demo

Single-page demo of `solveLambert` running in the browser via
`wasm-pack`'s `--target web` output.

## Build the wasm package

From the repo root:

```bash
wasm-pack build crates/lambert_izzo_wasm --target web --release \
    --out-dir examples/web/pkg
```

That populates `examples/web/pkg/` with the `.js`, `.d.ts`, and `.wasm`
files this demo imports.

## Serve

The browser will block ES-module imports from `file://` URLs, so serve
locally:

```bash
cd crates/lambert_izzo_wasm/examples/web
python3 -m http.server 8000
```

Then open <http://localhost:8000> and click **Solve**. Expected output:
the single-revolution `v1` / `v2` velocities and the Householder iteration
count from `response.diagnostics.single.iters`.

## What it shows

- Default request: a 90° circular LEO transfer at 7000 km altitude, Earth
  gravity (`mu = 398_600.4418 km³/s²`).
- The solve fires on the **Solve** button click — edit the form, then
  click to run `solveLambert` with the new request.
- Errors render the typed `error.kind` discriminator and message —
  trigger one by setting `r1` to `[0, 0, 0]` (degenerate) or `tof` to a
  negative number.
