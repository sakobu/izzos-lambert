//! Numerical tolerances and iteration limits for the Izzo Lambert solver.
//!
//! All literals used by the algorithm live here, no bare
//! `1e-X` should appear in the kernel modules.

// --- Geometry tolerances ---
//
// Guard near-degenerate inputs in the geometry construction stage.
// `MIN_POSITION_NORM_KM` is in km (matches the public API's `_km` convention).
// `COLINEARITY_TOL` is the unitless ratio `|r1 × r2| / (|r1| · |r2|)`.

/// Minimum norm for a valid position vector (km). Below this, the position is
/// treated as the origin and the chord/lambda construction is undefined.
///
/// `1e-15 km` (= 1 femtometer) is essentially `f64` noise on a unit-scale
/// vector. Practical Earth-orbit callers can rely on this catching only
/// truly-zero inputs from buggy upstream code.
pub const MIN_POSITION_NORM_KM: f64 = 1e-15;

/// Threshold below which `r1` and `r2` are treated as colinear. Compares
/// the unitless ratio `|r1 × r2| / (|r1| · |r2|)`, so the check has meaning
/// regardless of input scale.
pub const COLINEARITY_TOL: f64 = 1e-15;

// --- Convergence tolerances (Izzo §5) ---

/// Householder `|Δx|` tolerance for the single-revolution branch.
///
/// Izzo §5: "1e-5 is sufficient." Looser than the multi-rev tolerance
/// because the single-rev branch is well-conditioned; further refinement
/// does not improve the velocity reconstruction beyond `~1e-9` ΔE/E.
pub(crate) const HOUSEHOLDER_TOL_SINGLE: f64 = 1e-5;

/// Householder `|Δx|` tolerance for multi-revolution branches (`M ≥ 1`).
///
/// Tighter than single-rev because multi-rev branches lie close to the
/// `T_min` minimum where `dT/dx → 0` and the iteration is more sensitive.
pub(crate) const HOUSEHOLDER_TOL_MULTI: f64 = 1e-8;

/// Halley iteration tolerance for the `T_min` search on multi-rev branches.
///
/// Tighter than the Householder tolerances because `T_min` decides whether
/// a given `M` branch admits a solution at all — a noisy `T_min` could
/// silently drop or admit branches.
pub(crate) const HALLEY_TOL: f64 = 1e-13;

// --- Iteration limits ---

/// Householder safety cap. Izzo Fig. 3 reports convergence in `≤ 7` iters
/// over `1e6` random geometries; `15` is a comfortable margin.
pub(crate) const HOUSEHOLDER_MAX_ITERS: u32 = 15;

/// Halley safety cap. Cubic convergence reaches `f64` precision well under
/// this limit in normal use.
pub(crate) const HALLEY_MAX_ITERS: u32 = 12;

// --- Time-of-flight regime thresholds (Izzo Eq. 9, 18, 20) ---
//
// Selects the TOF formulation based on `|x − 1|`. The Lagrange form (Eq. 9)
// is numerically poor near the parabolic point `x = 1`; Lancaster (Eq. 18)
// is fine in a band; Battin's hypergeometric series (Eq. 20) is required
// closest to `x = 1`.

/// `|x − 1| ≤ this` → use Battin's hypergeometric formulation (Izzo Eq. 20).
///
/// Below `0.01`, the Lancaster form's `1/(1−x²)` factor amplifies round-off
/// to `~1e-6` of the true TOF; Battin's series remains exact.
pub(crate) const BATTIN_THRESHOLD: f64 = 0.01;

/// `BATTIN_THRESHOLD < |x − 1| ≤ this` → use Lancaster–Blanchard (Izzo Eq. 18).
///
/// Above `0.2`, the Lagrange form's trigonometric reductions are stable
/// enough that the Lancaster path's added complexity is unnecessary.
pub(crate) const LAGRANGE_THRESHOLD: f64 = 0.2;

// --- 2F1 hypergeometric series (Izzo Eq. 20) ---

/// Convergence tolerance for the direct `2F1(3, 1; 5/2; z)` series sum
/// used by the Battin TOF formulation.
///
/// `1e-11` is comfortably below the Battin regime's contribution to total
/// solver error and converges in `< 30` terms for `|z| < BATTIN_THRESHOLD`.
pub(crate) const HYPERGEOMETRIC_2F1_TOL: f64 = 1e-11;

/// Maximum terms for the `2F1` series sum before bailing out.
///
/// Reached only on pathological inputs (e.g. `NaN`); under normal conditions
/// the series converges at the tolerance above in `< 30` terms.
pub(crate) const HYPERGEOMETRIC_2F1_MAX_TERMS: u32 = 1000;
