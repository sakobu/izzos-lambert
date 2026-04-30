//! Read-only collection of multi-revolution solver diagnostics.
//!
//! Mirrors the structure of [`crate::MultiRevSet`] one-to-one: each entry's
//! `n_revs` and slot order matches the corresponding pair in `MultiRevSet`.
//! Populated alongside the solutions; consumers read it as a slice or
//! iterate.

use crate::MultiRevPairDiagnostics;
use crate::multi_rev_array::multi_rev_collection;

multi_rev_collection! {
    /// Per-pair diagnostics for the multi-revolution solver branches.
    ///
    /// Stack-allocated, capacity
    /// [`MAX_MULTI_REV_PAIRS`](crate::MAX_MULTI_REV_PAIRS), populated 1:1
    /// with [`crate::MultiRevSet`].
    MultiRevDiagnostics<MultiRevPairDiagnostics>;
    /// Number of entries.
    len();
    /// `true` if there are no diagnostics entries.
    is_empty();
    /// Iterate over the entries by reference, in ascending `M` order.
    iter();
}
