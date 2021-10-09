use std::fmt::Display;

use approx::{abs_diff_eq, AbsDiffEq};
use geo_types::{point, CoordNum, Point};
use num_traits::Float;

/// Represents the area outside of the triangulation.
/// Halfedges on the convex hull (which don't have an adjacent halfedge)
/// will have this value.
pub const EMPTY: usize = usize::max_value();

pub trait CoordType: CoordNum + AbsDiffEq<Epsilon = Self> + Display {
    fn infinity() -> Self;
    fn neg_infinity() -> Self;

    fn epsilon() -> Self;

    fn abs(self) -> Self;

    fn floor(self) -> Self;

    fn min(self, other: Self) -> Self;

    fn max(self, other: Self) -> Self;
}

pub trait DelaunayMath<T>
where
    T: CoordType,
{
    /// Near-duplicate points (where both `x` and `y` only differ within this value)
    /// will not be included in the triangulation for robustness.
    fn near_equals(a: Point<T>, b: Point<T>) -> bool {
        abs_diff_eq!(a, b, epsilon = T::epsilon())
    }

    fn circumdelta(a: Point<T>, b: Point<T>, c: Point<T>) -> Point<T> {
        let dx = b.x() - a.x();
        let dy = b.y() - a.y();
        let ex = c.x() - a.x();
        let ey = c.y() - a.y();

        let bl = dx * dx + dy * dy;
        let cl = ex * ex + ey * ey;
        let d = T::from(0.5).unwrap() / (dx * ey - dy * ex);

        let x = (ey * bl - dy * cl) * d;
        let y = (dx * cl - ex * bl) * d;

        point!(x: x, y: y)
    }

    fn in_circle(a: Point<T>, b: Point<T>, c: Point<T>, p: Point<T>) -> bool {
        let d = a - p;
        let e = b - p;
        let f = c - p;

        let ap = d.dot(d);
        let bp = e.dot(e);
        let cp = f.dot(f);

        d.x() * (e.y() * cp - bp * f.y()) - d.y() * (e.x() * cp - bp * f.x())
            + ap * (e.x() * f.y() - e.y() * f.x())
            < T::zero()
    }

    #[inline]
    fn circumradius2(a: Point<T>, b: Point<T>, c: Point<T>) -> T {
        let d = Self::circumdelta(a, b, c);
        d.dot(d)
    }

    #[inline]
    fn circumcenter(a: Point<T>, b: Point<T>, c: Point<T>) -> Point<T> {
        let d = Self::circumdelta(a, b, c);

        a + d
    }

    #[inline]
    fn dist2(p0: Point<T>, p: Point<T>) -> T {
        let d = p0 - p;
        d.dot(d)
    }

    #[inline]
    fn sortf(f: &mut Vec<(usize, T)>) {
        f.sort_unstable_by(|&(_, da), &(_, db)| da.partial_cmp(&db).unwrap());
    }

    #[inline]
    fn orient(p: Point<T>, q: Point<T>, r: Point<T>) -> bool {
        //(q.y() - p.y()) * (r.x() - q.x()) - (q.x() - p.x()) * (r.y() - q.y()) < 0.0
        p.cross_prod(q, r) >= T::zero()
    }

    fn find_closest_point(points: &[Point<T>], p0: Point<T>) -> Option<usize> {
        let mut min_dist = T::infinity();
        let mut k: usize = 0;
        for (i, p) in points.iter().enumerate() {
            let d = Self::dist2(p0, *p);
            if !d.is_zero() && d.lt(&min_dist) {
                k = i;
                min_dist = d;
            }
        }
        if min_dist == T::infinity() {
            None
        } else {
            Some(k)
        }
    }

    fn calc_bbox_center(points: &[Point<T>]) -> Point<T> {
        let mut min_x: T = T::infinity();
        let mut min_y: T = T::infinity();
        let mut max_x: T = T::neg_infinity();
        let mut max_y: T = T::neg_infinity();

        for p in points.iter() {
            min_x = min_x.min(p.x());
            min_y = min_y.min(p.y());
            max_x = max_x.max(p.x());
            max_y = max_y.max(p.y());
        }

        point!(
            x: (min_x + max_x) / T::from(2).unwrap(),
            y: (min_y + max_y) / T::from(2).unwrap()
        )
    }
}

impl CoordType for f32 {
    fn infinity() -> Self {
        Float::infinity()
    }

    fn neg_infinity() -> Self {
        Float::neg_infinity()
    }

    fn abs(self) -> Self {
        self.abs()
    }

    fn floor(self) -> Self {
        Float::abs(self)
    }

    fn min(self, other: Self) -> Self {
        self.min(other)
    }

    fn max(self, other: Self) -> Self {
        self.max(other)
    }

    fn epsilon() -> Self {
        f32::EPSILON * 2.
    }
}

impl CoordType for f64 {
    fn infinity() -> Self {
        Float::infinity()
    }

    fn neg_infinity() -> Self {
        Float::neg_infinity()
    }

    fn abs(self) -> Self {
        Float::abs(self)
    }

    fn floor(self) -> Self {
        self.floor()
    }

    fn min(self, other: Self) -> Self {
        self.min(other)
    }

    fn max(self, other: Self) -> Self {
        self.max(other)
    }

    fn epsilon() -> Self {
        f64::EPSILON * 2.
    }
}

impl CoordType for u64 {
    fn infinity() -> Self {
        u64::MAX
    }

    fn neg_infinity() -> Self {
        u64::MIN
    }

    fn floor(self) -> Self {
        self
    }

    fn min(self, other: Self) -> Self {
        Ord::min(self, other)
    }

    fn max(self, other: Self) -> Self {
        Ord::max(self, other)
    }

    fn abs(self) -> Self {
        self
    }

    fn epsilon() -> Self {
        2
    }
}

impl CoordType for usize {
    fn infinity() -> Self {
        usize::MAX
    }

    fn neg_infinity() -> Self {
        usize::MIN
    }

    fn abs(self) -> Self {
        self
    }

    fn floor(self) -> Self {
        self
    }

    fn min(self, other: Self) -> Self {
        Ord::min(self, other)
    }

    fn max(self, other: Self) -> Self {
        Ord::max(self, other)
    }

    fn epsilon() -> Self {
        2
    }
}

impl CoordType for i64 {
    fn infinity() -> Self {
        i64::MAX
    }

    fn neg_infinity() -> Self {
        i64::MIN
    }

    fn floor(self) -> Self {
        self
    }

    fn min(self, other: Self) -> Self {
        Ord::min(self, other)
    }

    fn max(self, other: Self) -> Self {
        Ord::max(self, other)
    }

    fn abs(self) -> Self {
        self
    }

    fn epsilon() -> Self {
        2
    }
}
