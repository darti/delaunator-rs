use std::fmt::Display;

use approx::{abs_diff_eq, AbsDiffEq};
use geo_types::{point, CoordFloat, CoordNum, Point};
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

    /// Near-duplicate points (where both `x` and `y` only differ within this value)
    /// will not be included in the triangulation for robustness.
    #[inline]
    fn near_equals(a: Point<Self>, b: Point<Self>) -> bool {
        abs_diff_eq!(a, b, epsilon = Self::epsilon())
    }

    #[inline]
    fn circumdelta(a: Point<Self>, b: Point<Self>, c: Point<Self>) -> Option<Point<Self>> {
        let d = b - a;
        let e = c - a;

        let bl = d.dot(d);
        let cl = e.dot(e);
        //let d = T::from(0.5).unwrap() / (dx * ey - dy * ex);

        let d_prime = Self::from(2).unwrap() * (d.x() * e.y() - d.y() * e.x());

        if d_prime.is_zero() {
            None
        } else {
            Some(
                point!(x: (e.y() * bl - d.y() * cl) / d_prime, y: (d.x() * cl - e.x() * bl) / d_prime),
            )
        }
    }

    #[inline]
    fn in_circle(a: Point<Self>, b: Point<Self>, c: Point<Self>, p: Point<Self>) -> bool {
        let d = a - p;
        let e = b - p;
        let f = c - p;

        let ap = d.dot(d);
        let bp = e.dot(e);
        let cp = f.dot(f);

        d.x() * (e.y() * cp - bp * f.y()) - d.y() * (e.x() * cp - bp * f.x())
            + ap * (e.x() * f.y() - e.y() * f.x())
            < Self::zero()
    }

    #[inline]
    fn circumradius2(a: Point<Self>, b: Point<Self>, c: Point<Self>) -> Option<Self> {
        Self::circumdelta(a, b, c).map(|d| d.dot(d))
    }

    #[inline]
    fn circumcenter(a: Point<Self>, b: Point<Self>, c: Point<Self>) -> Option<Point<Self>> {
        Self::circumdelta(a, b, c).map(|d| a + d)
    }

    #[inline]
    fn dist2(p0: Point<Self>, p: Point<Self>) -> Self {
        let d = p0 - p;
        d.dot(d)
    }

    #[inline]
    fn sortf(f: &mut Vec<(usize, Self)>) {
        f.sort_unstable_by(|&(_, da), &(_, db)| da.partial_cmp(&db).unwrap());
    }

    #[inline]
    fn orient(p: Point<Self>, q: Point<Self>, r: Point<Self>) -> bool {
        //(q.y() - p.y()) * (r.x() - q.x()) - (q.x() - p.x()) * (r.y() - q.y()) < 0.0
        p.cross_prod(q, r) >= Self::zero()
    }

    #[inline]
    fn find_closest_point(points: &[Point<Self>], p0: Point<Self>) -> Option<usize> {
        let mut min_dist = Self::infinity();
        let mut k: usize = 0;
        for (i, p) in points.iter().enumerate() {
            let d = Self::dist2(p0, *p);
            if !d.is_zero() && d.lt(&min_dist) {
                k = i;
                min_dist = d;
            }
        }
        if min_dist == Self::infinity() {
            None
        } else {
            Some(k)
        }
    }

    #[inline]
    fn calc_bbox_center(points: &[Point<Self>]) -> Point<Self> {
        let mut min_x = Self::infinity();
        let mut min_y = Self::infinity();
        let mut max_x = Self::neg_infinity();
        let mut max_y = Self::neg_infinity();

        for p in points.iter() {
            min_x = min_x.min(p.x());
            min_y = min_y.min(p.y());
            max_x = max_x.max(p.x());
            max_y = max_y.max(p.y());
        }

        point!(
            x: (min_x + max_x) / Self::from(2).unwrap(),
            y: (min_y + max_y) / Self::from(2).unwrap()
        )
    }
}

impl CoordType for f32 {
    #[inline]
    fn infinity() -> Self {
        Float::infinity()
    }

    #[inline]
    fn neg_infinity() -> Self {
        Float::neg_infinity()
    }

    #[inline]
    fn abs(self) -> Self {
        self.abs()
    }

    #[inline]
    fn floor(self) -> Self {
        Float::abs(self)
    }

    #[inline]
    fn min(self, other: Self) -> Self {
        self.min(other)
    }

    #[inline]
    fn max(self, other: Self) -> Self {
        self.max(other)
    }

    #[inline]
    fn epsilon() -> Self {
        f32::EPSILON * 2.
    }

    fn circumdelta(a: Point<Self>, b: Point<Self>, c: Point<Self>) -> Option<Point<Self>> {
        let d = b - a;
        let e = c - a;

        let bl = d.dot(d);
        let cl = e.dot(e);
        //let d = T::from(0.5).unwrap() / (dx * ey - dy * ex);

        let d_prime = 1. / (2. * (d.x() * e.y() - d.y() * e.x()));

        if d_prime == 0. {
            None
        } else {
            Some(
                point!(x: (e.y() * bl - d.y() * cl) * d_prime, y: (d.x() * cl - e.x() * bl) * d_prime),
            )
        }
    }
}

impl CoordType for f64 {
    #[inline]
    fn infinity() -> Self {
        Float::infinity()
    }

    #[inline]
    fn neg_infinity() -> Self {
        Float::neg_infinity()
    }

    #[inline]
    fn abs(self) -> Self {
        Float::abs(self)
    }

    #[inline]
    fn floor(self) -> Self {
        self.floor()
    }

    #[inline]
    fn min(self, other: Self) -> Self {
        self.min(other)
    }

    #[inline]
    fn max(self, other: Self) -> Self {
        self.max(other)
    }

    #[inline]
    fn epsilon() -> Self {
        f64::EPSILON * 2.
    }

    #[inline]
    fn circumdelta(a: Point<Self>, b: Point<Self>, c: Point<Self>) -> Option<Point<Self>> {
        let d = b - a;
        let e = c - a;

        let bl = d.dot(d);
        let cl = e.dot(e);
        //let d = T::from(0.5).unwrap() / (dx * ey - dy * ex);

        let d_prime = 1. / (2. * (d.x() * e.y() - d.y() * e.x()));

        if d_prime == 0. {
            None
        } else {
            Some(
                point!(x: (e.y() * bl - d.y() * cl) * d_prime, y: (d.x() * cl - e.x() * bl) * d_prime),
            )
        }
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
