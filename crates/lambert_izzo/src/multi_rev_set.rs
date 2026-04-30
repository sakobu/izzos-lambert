//! Read-only collection of multi-revolution solution pairs.
//!
//! Wraps an `ArrayVec<MultiRevPair, MAX_MULTI_REV_PAIRS>` so the underlying
//! `arrayvec` crate doesn't appear in the public API; consumers see only
//! `Deref<Target = [MultiRevPair]>` and `IntoIterator`. Capacity-bounded for
//! `no_std` use; population is crate-private.

use crate::MultiRevPair;
use crate::multi_rev_array::multi_rev_collection;

multi_rev_collection! {
    /// Multi-revolution Lambert solution pairs in ascending `M` order.
    ///
    /// Stack-allocated, capacity [`MAX_MULTI_REV_PAIRS`](crate::MAX_MULTI_REV_PAIRS).
    /// Populated by the solver; consumers read it as a slice
    /// (`Deref<Target = [MultiRevPair]>`) or iterate by value
    /// (`IntoIterator`). Empty for
    /// [`RevolutionBudget::SingleOnly`](crate::RevolutionBudget) and when
    /// no multi-rev branches are feasible for the given time of flight.
    MultiRevSet<MultiRevPair>;
    /// Number of multi-rev pairs in the set.
    len();
    /// `true` if no multi-rev pairs were found.
    is_empty();
    /// Iterate over the pairs by reference, in ascending `M` order.
    iter();
}
