use approx::AbsDiffEq;
use geo_types::{point, CoordFloat, Point};
use num_traits::{Float, NumCast};

pub trait CoordType: CoordFloat + AbsDiffEq<Epsilon = Self> {}

pub trait DelaunayMath<T>
where
    T: CoordType,
{
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

    #[inline(always)]
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
        let mut min_dist = Float::infinity();
        let mut k: usize = 0;

        for (i, p) in points.iter().enumerate() {
            let d = Self::dist2(p0, *p);
            if d > T::zero() && d < min_dist {
                k = i;
                min_dist = d;
            }
        }

        if min_dist == Float::infinity() {
            None
        } else {
            Some(k)
        }
    }

    fn calc_bbox_center(points: &[Point<T>]) -> Point<T> {
        let mut min_x: T = Float::infinity();
        let mut min_y: T = Float::infinity();
        let mut max_x: T = Float::neg_infinity();
        let mut max_y: T = Float::neg_infinity();

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

impl CoordType for f64 {}
impl CoordType for f32 {}
