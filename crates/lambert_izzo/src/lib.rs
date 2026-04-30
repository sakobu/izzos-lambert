//! Izzo's revisited Lambert solver — single + multi-revolution, short/long way.
//!
//! Reference: D. Izzo, *Revisiting Lambert's problem*, Celestial Mechanics &
//! Dynamical Astronomy, 2014. arXiv:1403.2705. PDF in `docs/izzo.pdf`.
//!
//! Inline `Eq. N` / `Algorithm N` references in the source point to that paper.
//!
//! # Public API
//!
//! Two entry points:
//!
//! - [`lambert`] — single Lambert solve from a [`LambertInput`].
//! - `lambert_par` (`rayon` feature) — parallel batch solve over a slice
//!   of [`LambertInput`]s.
//!
//! Both return [`LambertSolutions`], which always carries the
//! single-revolution trajectory, every reachable multi-rev pair, and the
//! per-branch [`SolverDiagnostics`] (Householder iteration counts).
//!
//! # Units
//!
//! The crate is unit-agnostic at the type level — all quantities are plain
//! `f64` and `[f64; 3]`. The convention used throughout the docs and
//! examples is SI for astrodynamics work:
//!
//! | Quantity                | Unit   |
//! | ----------------------- | ------ |
//! | Position                | km     |
//! | Velocity                | km/s   |
//! | Time of flight          | s      |
//! | Gravitational parameter | km³/s² |
//!
//! Any consistent unit system works — pass `r` in meters and `mu` in m³/s²
//! and you get velocities in m/s. The math is dimensionally homogeneous.
//!
//! The algorithm is also frame-invariant under any inertial frame — pass
//! `r1` and `r2` in the same inertial frame (ECI for Earth orbits, HCRS
//! for solar transfers, etc.) and the returned velocities are in that
//! same frame. The function signature is frame-agnostic; the calling
//! code's variable names carry the frame info.
//!
//! Position and velocity vectors are plain `[f64; 3]` arrays; the crate has
//! no hard math-library dependency. Both `nalgebra::Vector3<f64>` and
//! `glam::DVec3` already convert to/from `[f64; 3]` natively, so callers
//! using either library can pass and receive vectors without an explicit
//! interop layer.
//!
//! ```ignore
//! // nalgebra:
//! let r1: [f64; 3] = nalgebra::Vector3::new(7000.0, 0.0, 0.0).into();
//! let v1: nalgebra::Vector3<f64> = solution.single.v1.into();
//!
//! // glam:
//! let r2 = glam::DVec3::new(0.0, 7000.0, 0.0).to_array();
//! let v2 = glam::DVec3::from_array(solution.single.v2);
//! ```
//!
//! # Example
//!
//! ```
//! use lambert_izzo::{lambert, LambertError, LambertInput, RevolutionBudget, TransferWay};
//!
//! # fn main() -> Result<(), LambertError> {
//! // LEO → LEO 90° transfer at 7000 km altitude.
//! let mu = 398_600.4418;
//! let r = 7000.0_f64;
//! let input = LambertInput {
//!     r1: [r, 0.0, 0.0],
//!     r2: [0.0, r, 0.0],
//!     tof: core::f64::consts::PI / 2.0 * (r.powi(3) / mu).sqrt(),
//!     mu,
//!     way: TransferWay::Short,
//!     revolutions: RevolutionBudget::SingleOnly,
//! };
//!
//! let solutions = lambert(&input)?;
//! assert!(solutions.multi.is_empty());
//! assert!(solutions.diagnostics.single.iters > 0);
//! let v1 = solutions.single.v1;
//! # let _ = v1;
//! # Ok(())
//! # }
//! ```
//!
//! # Cargo features
//!
//! Both features are off by default.
//!
//! - **`serde`** — derives `Serialize` + `Deserialize` on every public
//!   type, including [`LambertError`]. Stays `no_std`-compatible.
//! - **`rayon`** — enables `lambert_par` for parallel batch evaluation.
//!   Pulls in `std` transitively and is **not** `no_std`-compatible.

#![cfg_attr(not(test), no_std)]
#![warn(missing_docs)]
#![warn(clippy::pedantic)]
#![warn(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::unreachable
)]
#![allow(clippy::module_name_repetitions)] // LambertError, LambertSolution, LambertSolutions, LambertDiagnostics

/// Hard upper bound on the number of multi-revolution pairs returned.
///
/// The Izzo formulation admits up to `⌊T/π⌋` multi-rev branches, which can
/// be arbitrarily large for very long times of flight — but practical
/// missions almost never exceed `M = 5`. This cap keeps the bounded return
/// type a fixed stack size.
///
/// The cap is type-enforced via [`BoundedRevs::MAX`]: requests above it
/// fail at construction time with [`RevsOutOfRange`] rather than being
/// silently truncated. The two constants are kept in sync by a static
/// assertion at the crate root.
pub const MAX_MULTI_REV_PAIRS: usize = 32;

// `BoundedRevs::MAX` (the type-level cap on requested revolutions) and
// `MAX_MULTI_REV_PAIRS` (the bounded-collection capacity holding the
// matching results) must agree so `RevolutionBudget::UpTo` is always
// representable in the return without truncation.
const _: () = assert!(BoundedRevs::MAX as usize == MAX_MULTI_REV_PAIRS); // u32 → usize: always safe (usize ≥ 32 bits)

/// Direction around the transfer plane from `r1` to `r2`.
///
/// `Short` is the geodesic transfer (`θ ≤ π`); `Long` traverses the other
/// way (`θ > π`). This is independent of prograde/retrograde — the orbit's
/// angular-momentum direction is set by the order of the `(r1, r2)` arguments.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum TransferWay {
    /// `θ ≤ π` — the short geodesic arc.
    Short,
    /// `θ > π` — the long way around the transfer plane.
    Long,
}

/// Maximum number of complete revolutions to consider beyond single-rev.
///
/// Multi-revolution branches admit two solutions per revolution count `M`
/// (long-period and short-period), so the total solution count is
/// `1 + 2 · min(max(), ⌊T/π⌋)` adjusted downward when a branch's `T_min`
/// exceeds the requested time of flight. The upper bound on `M` is
/// type-enforced via [`BoundedRevs`] — out-of-range requests fail at
/// construction time rather than being clamped.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum RevolutionBudget {
    /// Solve the single-revolution case only — always exactly one solution.
    #[default]
    SingleOnly,
    /// Search up to `M` complete revolutions inclusive, with `M` validated
    /// at construction time to lie in `1..=BoundedRevs::MAX`.
    UpTo(BoundedRevs),
}

impl RevolutionBudget {
    /// Wrap a pre-validated [`BoundedRevs`] into a budget. Infallible — the
    /// `1..=BoundedRevs::MAX` bound is already established at the
    /// [`BoundedRevs::try_new`] call site.
    #[must_use]
    pub const fn up_to(revs: BoundedRevs) -> Self {
        Self::UpTo(revs)
    }

    /// Ergonomic fallible constructor.
    ///
    /// Equivalent to `BoundedRevs::try_new(n).map(RevolutionBudget::up_to)`.
    /// Use [`RevolutionBudget::SingleOnly`] directly when you want to skip
    /// multi-rev entirely — passing `0` is rejected.
    ///
    /// # Errors
    ///
    /// Returns [`RevsOutOfRange`] when `n == 0` or `n > BoundedRevs::MAX`.
    pub const fn try_up_to(n: u32) -> Result<Self, RevsOutOfRange> {
        match BoundedRevs::try_new(n) {
            Ok(revs) => Ok(Self::UpTo(revs)),
            Err(e) => Err(e),
        }
    }

    /// The maximum revolution count this budget will search.
    ///
    /// Returns [`None`] for [`RevolutionBudget::SingleOnly`] and
    /// `Some(b)` for [`RevolutionBudget::UpTo`]. For driving a multi-rev
    /// loop, prefer [`RevolutionBudget::iter_revs`] — `max()` is intended
    /// for observability (display, logging, capacity hints).
    #[must_use]
    pub const fn max(self) -> Option<BoundedRevs> {
        match self {
            Self::SingleOnly => None,
            Self::UpTo(b) => Some(b),
        }
    }

    /// Iterator over the multi-rev branch counts this budget will search.
    ///
    /// Yields validated [`BoundedRevs`] values `1..=b` for
    /// [`RevolutionBudget::UpTo(b)`](RevolutionBudget::UpTo) and is empty
    /// for [`RevolutionBudget::SingleOnly`]. This is the canonical way to
    /// drive the multi-rev loop in the kernel — emitting `BoundedRevs`
    /// rather than raw `u32` keeps the `1..=BoundedRevs::MAX` invariant
    /// inside the type system all the way to the
    /// [`MultiRevPair::n_revs`] / [`MultiRevPairDiagnostics::n_revs`]
    /// fields on the returned [`LambertSolutions`]. Call `.get()` at any
    /// site that needs the raw count for arithmetic.
    pub fn iter_revs(self) -> impl Iterator<Item = BoundedRevs> {
        self.max()
            .into_iter()
            .flat_map(BoundedRevs::range_inclusive_one_to_self)
    }
}

mod bounded_revs;
mod constants;
mod error;
mod geometry;
mod multi_rev;
mod root_finding;
mod tof;
mod vec3;

#[cfg(test)]
mod tests;

pub use bounded_revs::{BoundedRevs, RevsOutOfRange};
pub use error::{LambertError, NonFiniteParameter, Position};
pub use multi_rev::{MultiRevDiagnostics, MultiRevSet};

use geometry::Geometry;
use root_finding::{Root, find_xy};

/// One Lambert transfer trajectory.
///
/// Pure trajectory data — start and end velocities. Solver diagnostics
/// live alongside the trajectory in [`LambertSolutions::diagnostics`].
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct LambertSolution {
    /// Velocity at `r1` (same units and frame as the inputs).
    pub v1: [f64; 3],
    /// Velocity at `r2` (same units and frame as the inputs).
    pub v2: [f64; 3],
}

/// One multi-revolution pair: long-period and short-period trajectories
/// for a given revolution count.
///
/// The Izzo formulation admits exactly two trajectories for each `M ≥ 1`:
/// the long-period branch (more time near apoapsis) and the short-period
/// branch (more time near periapsis).
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct MultiRevPair {
    /// Branch revolution count.
    pub n_revs: BoundedRevs,
    /// Long-period trajectory.
    pub long_period: LambertSolution,
    /// Short-period trajectory.
    pub short_period: LambertSolution,
}

/// All Lambert trajectories for a given boundary problem and revolution budget.
///
/// Always carries:
/// - `single` — the single-revolution trajectory;
/// - `multi` — every reachable multi-rev pair in ascending `M` order,
///   bounded by [`MAX_MULTI_REV_PAIRS`];
/// - `diagnostics` — per-branch [`SolverDiagnostics`] (Householder iteration
///   counts), structurally parallel to `single` / `multi`.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct LambertSolutions {
    /// Single-revolution trajectory — always present.
    pub single: LambertSolution,
    /// Multi-revolution pairs in ascending `M` order; empty for
    /// [`RevolutionBudget::SingleOnly`] or when no multi-rev branches are
    /// feasible for the given time of flight.
    pub multi: MultiRevSet,
    /// Per-branch solver diagnostics, structurally parallel to
    /// `single` / `multi`.
    pub diagnostics: LambertDiagnostics,
}

impl LambertSolutions {
    /// Highest revolution count for which the solver found a feasible
    /// `(long_period, short_period)` pair at the requested time-of-flight
    /// — the silent-skip boundary referenced in [`lambert`]'s
    /// _Validity / near-degenerate behavior_ section (the paper's `M`).
    ///
    /// Returns [`None`] for [`RevolutionBudget::SingleOnly`] or when no
    /// multi-rev branch was feasible at the requested TOF; otherwise
    /// `Some(b)` equal to `self.multi.last().unwrap().n_revs`. Pairs
    /// structurally with [`RevolutionBudget::max`] — both expose "highest
    /// relevant `M`" through the same [`BoundedRevs`] type.
    #[must_use]
    pub fn max_feasible_revs(&self) -> Option<BoundedRevs> {
        self.multi.last().map(|p| p.n_revs)
    }
}

/// Diagnostic data for one converged Householder solve.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct SolverDiagnostics {
    /// Householder iterations used to converge.
    pub iters: u32,
}

/// Diagnostics for one multi-rev pair.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct MultiRevPairDiagnostics {
    /// Branch revolution count.
    pub n_revs: BoundedRevs,
    /// Long-period branch diagnostics.
    pub long_period: SolverDiagnostics,
    /// Short-period branch diagnostics.
    pub short_period: SolverDiagnostics,
}

/// Diagnostics structure mirroring [`LambertSolutions`].
///
/// Carried inline by every [`LambertSolutions`]; the `single` / `multi`
/// layout matches the corresponding solution fields one-to-one.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct LambertDiagnostics {
    /// Single-rev solver diagnostics.
    pub single: SolverDiagnostics,
    /// Multi-rev pair diagnostics in the same order as
    /// [`LambertSolutions::multi`].
    pub multi: MultiRevDiagnostics,
}

/// One Lambert call's inputs.
///
/// All scalar fields use the crate's documented unit convention (km, s,
/// km/s, km³/s²) when read in SI; any consistent unit system works.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct LambertInput {
    /// Initial position, any consistent inertial frame.
    pub r1: [f64; 3],
    /// Final position, same frame as `r1`.
    pub r2: [f64; 3],
    /// Time of flight, `> 0`.
    pub tof: f64,
    /// Gravitational parameter of the central body, `> 0`.
    pub mu: f64,
    /// Short or long way around the transfer plane.
    pub way: TransferWay,
    /// Revolution budget — see [`RevolutionBudget`].
    pub revolutions: RevolutionBudget,
}

/// Solve Lambert's boundary-value problem using Izzo's revisited algorithm.
///
/// Householder iteration over Lancaster's free parameter `x`, dispatching
/// across three TOF regimes (Battin / Lancaster–Blanchard / Lagrange) for
/// numerical stability. Mathematically frame-invariant — pass `r1` and
/// `r2` in any consistent inertial frame and the returned velocities are
/// in that same frame.
///
/// Returns the always-present single-revolution trajectory, every
/// reachable multi-rev branch up to the configured [`RevolutionBudget`]
/// cap, and the per-branch
/// [`SolverDiagnostics`]. The revolution budget is itself capped at
/// [`BoundedRevs::MAX`] at construction time, so the returned `multi` set
/// fits within [`MAX_MULTI_REV_PAIRS`] without runtime clamping.
///
/// All quantities are dimensionally homogeneous — pass any consistent
/// unit system. The crate's docs and examples use km / km/s / s / km³/s².
///
/// # Invariants
///
/// All preconditions are validated at entry and returned as `Err(...)` on
/// violation — never panicked.
///
/// - `input.tof > 0`
/// - `input.mu > 0`
/// - `input.r1`, `input.r2`, `input.tof`, and `input.mu` are finite.
/// - `|input.r1|` and `|input.r2|` exceed an internal floor (`1e-15`).
/// - Transfer angle ∉ {0, π}, equivalently
///   `|r1 × r2| / (|r1| · |r2|)` exceeds an internal floor (`1e-15`).
///
/// # Validity / near-degenerate behavior
///
/// - **Transfer angle near `0` or `π`** — the transfer plane is undefined;
///   returns [`LambertError::CollinearGeometry`]. Callers near these
///   boundaries should perturb one position by a tiny off-plane offset.
/// - **Near-parabolic** — the Lagrange and Lancaster TOF formulations lose
///   precision near `x = 1`; the solver switches to Battin's hypergeometric
///   series (Izzo Eq. 20) automatically.
/// - **Hyperbolic transfers (`x > 1`)** — admitted on the single-rev branch;
///   multi-revolution solutions do not exist on a hyperbola and are silently
///   skipped.
/// - **Multi-rev infeasibility** — for `M ≥ 1`, the branch admits a solution
///   only when `tof ≥ T_min(M, λ)`. Higher-`M` branches are dropped from the
///   returned `multi` set when their `T_min` exceeds the requested TOF;
///   call [`LambertSolutions::max_feasible_revs`] to detect this boundary
///   programmatically.
///
/// # Returns
///
/// [`LambertSolutions`] with the always-present `single` trajectory, a
/// `multi` set of [`MultiRevPair`] entries in ascending `M` order, and the
/// matching [`LambertDiagnostics`]. The `multi` set is empty for
/// [`RevolutionBudget::SingleOnly`] and may be shorter than the configured
/// cap if higher-`M` branches are infeasible — in which case
/// [`LambertSolutions::max_feasible_revs`] returns the highest `M` actually
/// solved. Each [`LambertSolution`] carries the start and end velocities
/// `(v1, v2)` in the input frame and units.
///
/// # Errors
///
/// - [`LambertError::NonFiniteInput`] — any scalar input or position vector
///   component is `NaN`, `+inf`, or `-inf`.
/// - [`LambertError::NonPositiveTimeOfFlight`] — `tof <= 0`.
/// - [`LambertError::NonPositiveMu`] — `mu <= 0`.
/// - [`LambertError::DegeneratePositionVector`] — `|r1|` or `|r2|` below the
///   internal position-norm floor.
/// - [`LambertError::CollinearGeometry`] — transfer plane undefined.
/// - [`LambertError::NoConvergence`] / [`LambertError::SingularDenominator`]
///   — Householder iteration failed.
///
/// # Examples
///
/// Multi-revolution search up to `M = 3` on a 90° LEO transfer with a
/// time of flight long enough for at least one revolution to be feasible:
///
/// ```
/// use lambert_izzo::{lambert, LambertInput, RevolutionBudget, TransferWay};
///
/// # fn main() -> Result<(), Box<dyn std::error::Error>> {
/// let mu = 398_600.4418;
/// let r = 7000.0_f64;
/// let period = 2.0 * core::f64::consts::PI * (r.powi(3) / mu).sqrt();
///
/// let solutions = lambert(&LambertInput {
///     r1: [r, 0.0, 0.0],
///     r2: [0.0, r, 0.0],
///     tof: 2.5 * period,                              // long enough for M ≥ 1
///     mu,
///     way: TransferWay::Short,
///     revolutions: RevolutionBudget::try_up_to(3)?,   // RevsOutOfRange if > MAX
/// })?;
///
/// assert!(!solutions.multi.is_empty());
/// assert!(solutions.multi.len() <= 3);
/// # Ok(())
/// # }
/// ```
pub fn lambert(input: &LambertInput) -> Result<LambertSolutions, LambertError> {
    let geom = Geometry::from_inputs(input.r1, input.r2, input.tof, input.mu, input.way)?;
    let roots = find_xy(&geom, input.revolutions)?;
    Ok(build_solutions(&geom, &roots))
}

/// Parallel batch solver over a slice of [`LambertInput`]s (Rayon-backed).
///
/// Yields one `Result<LambertSolutions, LambertError>` per input. Output
/// ordering is preserved when the caller consumes via `.collect::<Vec<_>>()`
/// or any indexed Rayon adaptor; reduce-style consumers (`.for_each`,
/// `.sum`) see results in arbitrary order.
///
/// # Invariants
///
/// Per-input — same as [`lambert`].
///
/// # Validity / near-degenerate behavior
///
/// Per-input — same as [`lambert`].
///
/// # Errors
///
/// Each yielded item is independent — one input failing does not stop the
/// other parallel tasks. Per-input error variants are the same as [`lambert`].
#[cfg(feature = "rayon")]
#[allow(clippy::must_use_candidate)] // Rayon's IndexedParallelIterator isn't #[must_use]; the caller is expected to chain `.for_each` / `.collect`.
pub fn lambert_par(
    inputs: &[LambertInput],
) -> impl rayon::iter::IndexedParallelIterator<Item = Result<LambertSolutions, LambertError>> + '_ {
    use rayon::prelude::*;
    inputs.par_iter().map(lambert)
}

#[allow(clippy::similar_names)] // v_r1/v_r2/v_t1/v_t2 are radial/tangential velocity components at points 1/2 — Izzo §2.
fn reconstruct(geom: &Geometry, root: &Root) -> LambertSolution {
    // Velocity reconstruction (Izzo Algorithm 1).
    // Paper-named locals; `lambda*y ± x` are the elliptic-anomaly combos
    // that appear repeatedly in the radial/tangential decomposition.
    let lambda_y_minus_x = geom.lambda * root.y - root.x;
    let lambda_y_plus_x = geom.lambda * root.y + root.x;
    let tangential_num = geom.gamma * geom.sigma * (root.y + geom.lambda * root.x);

    let v_r1 = geom.gamma * (lambda_y_minus_x - geom.rho * lambda_y_plus_x) / geom.r1n;
    let v_r2 = -geom.gamma * (lambda_y_minus_x + geom.rho * lambda_y_plus_x) / geom.r2n;
    let v_t1 = tangential_num / geom.r1n;
    let v_t2 = tangential_num / geom.r2n;

    // v_radial · ir + v_tangential · it, component-wise.
    let v1 = [
        geom.ir1[0] * v_r1 + geom.it1[0] * v_t1,
        geom.ir1[1] * v_r1 + geom.it1[1] * v_t1,
        geom.ir1[2] * v_r1 + geom.it1[2] * v_t1,
    ];
    let v2 = [
        geom.ir2[0] * v_r2 + geom.it2[0] * v_t2,
        geom.ir2[1] * v_r2 + geom.it2[1] * v_t2,
        geom.ir2[2] * v_r2 + geom.it2[2] * v_t2,
    ];
    LambertSolution { v1, v2 }
}

fn build_solutions(geom: &Geometry, roots: &root_finding::Roots) -> LambertSolutions {
    let single = reconstruct(geom, &roots.single);
    let single_diagnostics = SolverDiagnostics {
        iters: roots.single.iters,
    };

    let mut multi = MultiRevSet::new();
    let mut multi_diagnostics = MultiRevDiagnostics::new();
    for pair in &roots.multi {
        // `try_push` (not `push`) because lib code is panic-free; `roots.multi`
        // is itself bounded at `MAX_MULTI_REV_PAIRS`, so the push cannot fail.
        multi.try_push(MultiRevPair {
            n_revs: pair.n_revs,
            long_period: reconstruct(geom, &pair.long_period),
            short_period: reconstruct(geom, &pair.short_period),
        });
        multi_diagnostics.try_push(MultiRevPairDiagnostics {
            n_revs: pair.n_revs,
            long_period: SolverDiagnostics {
                iters: pair.long_period.iters,
            },
            short_period: SolverDiagnostics {
                iters: pair.short_period.iters,
            },
        });
    }
    LambertSolutions {
        single,
        multi,
        diagnostics: LambertDiagnostics {
            single: single_diagnostics,
            multi: multi_diagnostics,
        },
    }
}
