//! Read-only collections of multi-revolution solver outputs.
//!
//! [`MultiRevSet`] holds the solution pairs and [`MultiRevDiagnostics`] holds
//! the per-pair convergence stats; both wrap an
//! `ArrayVec<_, MAX_MULTI_REV_PAIRS>` so the underlying `arrayvec` crate
//! doesn't appear in the public API. Consumers see only
//! `Deref<Target = [_]>` and `IntoIterator`. Capacity-bounded for `no_std`
//! use; population is crate-private.
//!
//! The two types share an identical surface (`len`, `is_empty`, `iter`,
//! deref-to-slice, by-value and by-ref `IntoIterator`) — the
//! `multi_rev_collection!` macro defined below generates both so each can
//! carry its own rustdoc and namespacing without duplicating the impls.

use crate::{MultiRevPair, MultiRevPairDiagnostics};

macro_rules! multi_rev_collection {
    (
        $(#[$type_doc:meta])*
        $name:ident < $elem:ty >;
        $(#[$len_doc:meta])*
        len();
        $(#[$is_empty_doc:meta])*
        is_empty();
        $(#[$iter_doc:meta])*
        iter();
    ) => {
        $(#[$type_doc])*
        #[derive(Debug, Clone, Default, PartialEq)]
        #[cfg_attr(
            feature = "serde",
            derive(serde::Serialize, serde::Deserialize),
            serde(transparent)
        )]
        pub struct $name {
            inner: arrayvec::ArrayVec<$elem, { $crate::MAX_MULTI_REV_PAIRS }>,
        }

        impl $name {
            pub(crate) fn new() -> Self {
                Self {
                    inner: arrayvec::ArrayVec::new(),
                }
            }

            pub(crate) fn try_push(&mut self, value: $elem) -> bool {
                self.inner.try_push(value).is_ok()
            }

            $(#[$len_doc])*
            #[must_use]
            pub fn len(&self) -> usize {
                self.inner.len()
            }

            $(#[$is_empty_doc])*
            #[must_use]
            pub fn is_empty(&self) -> bool {
                self.inner.is_empty()
            }

            $(#[$iter_doc])*
            pub fn iter(&self) -> core::slice::Iter<'_, $elem> {
                self.inner.iter()
            }
        }

        impl core::ops::Deref for $name {
            type Target = [$elem];
            fn deref(&self) -> &[$elem] {
                &self.inner
            }
        }

        impl IntoIterator for $name {
            type Item = $elem;
            type IntoIter =
                arrayvec::IntoIter<$elem, { $crate::MAX_MULTI_REV_PAIRS }>;
            fn into_iter(self) -> Self::IntoIter {
                self.inner.into_iter()
            }
        }

        impl<'a> IntoIterator for &'a $name {
            type Item = &'a $elem;
            type IntoIter = core::slice::Iter<'a, $elem>;
            fn into_iter(self) -> Self::IntoIter {
                self.inner.iter()
            }
        }
    };
}

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
