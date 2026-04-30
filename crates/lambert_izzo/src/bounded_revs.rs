//! Type-enforced revolution count for [`crate::RevolutionBudget::UpTo`].
//!
//! [`BoundedRevs`] is a newtype around [`NonZeroU32`] that additionally
//! constrains the value to `1..=BoundedRevs::MAX`. The cap mirrors
//! [`crate::MAX_MULTI_REV_PAIRS`] (the bounded-collection capacity used by
//! the solver's return types) so that any constructed budget is, by
//! construction, representable in the output without truncation.
//!
//! Construction is fallible — out-of-range requests fail loudly via
//! [`RevsOutOfRange`] rather than being silently clamped.

use core::num::NonZeroU32;

/// Validated revolution count for multi-rev Lambert solves.
///
/// Constructed via [`BoundedRevs::try_new`], which enforces the
/// `1..=BoundedRevs::MAX` invariant. Use [`crate::RevolutionBudget::up_to`]
/// (total) or [`crate::RevolutionBudget::try_up_to`] (fallible) to wrap the
/// validated value into the budget enum consumed by [`crate::lambert`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct BoundedRevs(NonZeroU32);

impl BoundedRevs {
    /// Hard upper bound on the revolution count.
    ///
    /// Mirrors [`crate::MAX_MULTI_REV_PAIRS`]; the link is static-asserted
    /// at the crate root.
    pub const MAX: u32 = 32;

    /// Construct a [`BoundedRevs`] from a raw `u32`.
    ///
    /// # Errors
    ///
    /// Returns [`RevsOutOfRange`] when `n == 0` or `n > BoundedRevs::MAX`.
    pub const fn try_new(n: u32) -> Result<Self, RevsOutOfRange> {
        if let Some(nz) = NonZeroU32::new(n) {
            if nz.get() <= Self::MAX {
                return Ok(Self(nz));
            }
        }
        Err(RevsOutOfRange {
            requested: n,
            max: Self::MAX,
        })
    }

    /// The wrapped revolution count (`1..=BoundedRevs::MAX`).
    #[must_use]
    pub const fn get(self) -> u32 {
        self.0.get()
    }

    /// Iterator over `1..=self` as validated [`BoundedRevs`] values.
    ///
    /// Each yielded value inherits the `1..=BoundedRevs::MAX` invariant from
    /// `self`, so the bounds check inside [`BoundedRevs::try_new`] would be
    /// redundant — values are constructed via the private tuple constructor
    /// instead. Used by [`crate::RevolutionBudget::iter_revs`].
    pub(crate) fn range_inclusive_one_to_self(self) -> impl Iterator<Item = Self> {
        (1..=self.0.get()).filter_map(|n| NonZeroU32::new(n).map(Self))
    }
}

impl core::fmt::Display for BoundedRevs {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        self.0.fmt(f)
    }
}

#[cfg(feature = "serde")]
impl serde::Serialize for BoundedRevs {
    fn serialize<S: serde::Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        self.0.get().serialize(serializer)
    }
}

#[cfg(feature = "serde")]
impl<'de> serde::Deserialize<'de> for BoundedRevs {
    fn deserialize<D: serde::Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        let n = u32::deserialize(deserializer)?;
        Self::try_new(n).map_err(serde::de::Error::custom)
    }
}

/// Construction-time validation error returned by [`BoundedRevs::try_new`]
/// (and by extension [`crate::RevolutionBudget::try_up_to`]).
///
/// Distinct from [`crate::LambertError`] — that type carries solver-runtime
/// failures; this one carries the precondition violation that prevents the
/// solver from being called in the first place.
#[derive(Debug, Clone, Copy, PartialEq, Eq, thiserror::Error)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[error("revolution count {requested} out of range (must be 1..={max})")]
pub struct RevsOutOfRange {
    /// The rejected value the caller passed.
    pub requested: u32,
    /// The inclusive upper bound (`BoundedRevs::MAX`).
    pub max: u32,
}

#[cfg(test)]
mod tests {
    #![allow(clippy::unwrap_used)]

    use super::*;

    #[test]
    fn try_new_accepts_one_through_max() {
        for n in 1..=BoundedRevs::MAX {
            let revs = BoundedRevs::try_new(n).unwrap();
            assert_eq!(revs.get(), n);
        }
    }

    #[test]
    fn try_new_rejects_zero() {
        let err = BoundedRevs::try_new(0).unwrap_err();
        assert_eq!(
            err,
            RevsOutOfRange {
                requested: 0,
                max: BoundedRevs::MAX,
            }
        );
    }

    #[test]
    fn try_new_rejects_max_plus_one() {
        let err = BoundedRevs::try_new(BoundedRevs::MAX + 1).unwrap_err();
        assert_eq!(
            err,
            RevsOutOfRange {
                requested: BoundedRevs::MAX + 1,
                max: BoundedRevs::MAX,
            }
        );
    }

    #[test]
    fn try_new_rejects_u32_max() {
        let err = BoundedRevs::try_new(u32::MAX).unwrap_err();
        assert_eq!(
            err,
            RevsOutOfRange {
                requested: u32::MAX,
                max: BoundedRevs::MAX,
            }
        );
    }

    #[test]
    fn try_new_is_const() {
        // Compiles only because `try_new` is `const fn` — regression guard.
        const REVS: Result<BoundedRevs, RevsOutOfRange> = BoundedRevs::try_new(5);
        assert_eq!(REVS.unwrap().get(), 5);
    }

    #[test]
    fn display_matches_inner_value() {
        let revs = BoundedRevs::try_new(5).unwrap();
        assert_eq!(format!("{revs}"), "5");
    }

    #[cfg(feature = "serde")]
    #[test]
    fn serde_round_trip_preserves_value() {
        let revs = BoundedRevs::try_new(7).unwrap();
        let json = serde_json::to_string(&revs).unwrap();
        assert_eq!(json, "7");
        let back: BoundedRevs = serde_json::from_str(&json).unwrap();
        assert_eq!(revs, back);
    }

    #[cfg(feature = "serde")]
    #[test]
    fn serde_rejects_out_of_range_on_deserialize() {
        // 0 and > MAX must not slip through deserialization.
        assert!(serde_json::from_str::<BoundedRevs>("0").is_err());
        let too_big = format!("{}", BoundedRevs::MAX + 1);
        assert!(serde_json::from_str::<BoundedRevs>(&too_big).is_err());
    }
}
