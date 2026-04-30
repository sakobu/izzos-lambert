//! Internal macro behind [`MultiRevSet`](crate::MultiRevSet) and
//! [`MultiRevDiagnostics`](crate::MultiRevDiagnostics). Both expose the
//! same read-only, capacity-bounded `ArrayVec` view; the macro lets each
//! type carry its own rustdoc and namespacing while sharing a single
//! implementation.

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

pub(crate) use multi_rev_collection;
