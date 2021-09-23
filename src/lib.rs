/*!
A very fast 2D [Delaunay Triangulation](https://en.wikipedia.org/wiki/Delaunay_triangulation) library for Rust.
A port of [Delaunator](https://github.com/mapbox/delaunator).

# Example

```rust
use delaunator::{Point, triangulate};

let points = vec![
    Point { x: 0., y: 0. },
    Point { x: 1., y: 0. },
    Point { x: 1., y: 1. },
    Point { x: 0., y: 1. },
];

let result = triangulate(&points);
println!("{:?}", result.triangles); // [0, 2, 1, 0, 3, 2]
```
*/

use approx::{abs_diff_eq, AbsDiffEq};
use num_traits::{Float, NumCast};

use geo_types::{point, CoordFloat, MultiPoint, Point};

/// Near-duplicate points (where both `x` and `y` only differ within this value)
/// will not be included in the triangulation for robustness.
pub const EPSILON: f64 = f64::EPSILON * 2.0;

fn dist2<T>(p0: Point<T>, p: &Point<T>) -> T
where
    T: CoordFloat,
{
    let d = p0 - *p;
    d.dot(d)
}

fn orient<T>(p: Point<T>, q: Point<T>, r: Point<T>) -> bool
where
    T: CoordFloat,
{
    //(q.y() - p.y()) * (r.x() - q.x()) - (q.x() - p.x()) * (r.y() - q.y()) < 0.0
    p.cross_prod(q, r) >= T::zero()
}

fn circumdelta<T>(a: Point<T>, b: Point<T>, c: Point<T>) -> Point<T>
where
    T: CoordFloat,
{
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

fn circumradius2<T>(a: Point<T>, b: Point<T>, c: Point<T>) -> T
where
    T: CoordFloat,
{
    let d = circumdelta(a, b, c);
    d.dot(d)
}

fn circumcenter<T>(a: Point<T>, b: Point<T>, c: Point<T>) -> Point<T>
where
    T: CoordFloat,
{
    let d = circumdelta(a, b, c);

    a + d
}

fn in_circle<T>(a: Point<T>, b: Point<T>, c: Point<T>, p: Point<T>) -> bool
where
    T: CoordFloat,
{
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

/// Represents the area outside of the triangulation.
/// Halfedges on the convex hull (which don't have an adjacent halfedge)
/// will have this value.
pub const EMPTY: usize = usize::max_value();

/// Next halfedge in a triangle.
pub fn next_halfedge(i: usize) -> usize {
    if i % 3 == 2 {
        i - 2
    } else {
        i + 1
    }
}

/// Previous halfedge in a triangle.
pub fn prev_halfedge(i: usize) -> usize {
    if i % 3 == 0 {
        i + 2
    } else {
        i - 1
    }
}

/// Result of the Delaunay triangulation.
pub struct Triangulation {
    /// A vector of point indices where each triple represents a Delaunay triangle.
    /// All triangles are directed counter-clockwise.
    pub triangles: Vec<usize>,

    /// A vector of adjacent halfedge indices that allows traversing the triangulation graph.
    ///
    /// `i`-th half-edge in the array corresponds to vertex `triangles[i]`
    /// the half-edge is coming from. `halfedges[i]` is the index of a twin half-edge
    /// in an adjacent triangle (or `EMPTY` for outer half-edges on the convex hull).
    pub halfedges: Vec<usize>,

    /// A vector of indices that reference points on the convex hull of the triangulation,
    /// counter-clockwise.
    pub hull: Vec<usize>,
}

impl Triangulation {
    fn new(n: usize) -> Self {
        let max_triangles = if n > 2 { 2 * n - 5 } else { 0 };

        Self {
            triangles: Vec::with_capacity(max_triangles * 3),
            halfedges: Vec::with_capacity(max_triangles * 3),
            hull: Vec::new(),
        }
    }

    /// The number of triangles in the triangulation.
    pub fn len(&self) -> usize {
        self.triangles.len() / 3
    }

    fn add_triangle(
        &mut self,
        i0: usize,
        i1: usize,
        i2: usize,
        a: usize,
        b: usize,
        c: usize,
    ) -> usize {
        let t = self.triangles.len();

        self.triangles.push(i0);
        self.triangles.push(i1);
        self.triangles.push(i2);

        self.halfedges.push(a);
        self.halfedges.push(b);
        self.halfedges.push(c);

        if a != EMPTY {
            self.halfedges[a] = t;
        }
        if b != EMPTY {
            self.halfedges[b] = t + 1;
        }
        if c != EMPTY {
            self.halfedges[c] = t + 2;
        }

        t
    }

    fn legalize<T>(&mut self, a: usize, points: &[Point<T>], hull: &mut Hull<T>) -> usize
    where
        T: CoordFloat,
    {
        let b = self.halfedges[a];

        // if the pair of triangles doesn't satisfy the Delaunay condition
        // (p1 is inside the circumcircle of [p0, pl, pr]), flip them,
        // then do the same check/flip recursively for the new pair of triangles
        //
        //           pl                    pl
        //          /||\                  /  \
        //       al/ || \bl            al/    \a
        //        /  ||  \              /      \
        //       /  a||b  \    flip    /___ar___\
        //     p0\   ||   /p1   =>   p0\---bl---/p1
        //        \  ||  /              \      /
        //       ar\ || /br             b\    /br
        //          \||/                  \  /
        //           pr                    pr
        //
        let ar = prev_halfedge(a);

        if b == EMPTY {
            return ar;
        }

        let al = next_halfedge(a);
        let bl = prev_halfedge(b);

        let p0 = self.triangles[ar];
        let pr = self.triangles[a];
        let pl = self.triangles[al];
        let p1 = self.triangles[bl];

        let illegal = in_circle(points[p0], points[pr], points[pl], points[p1]);
        if illegal {
            self.triangles[a] = p1;
            self.triangles[b] = p0;

            let hbl = self.halfedges[bl];
            let har = self.halfedges[ar];

            // edge swapped on the other side of the hull (rare); fix the halfedge reference
            if hbl == EMPTY {
                let mut e = hull.start;
                loop {
                    if hull.tri[e] == bl {
                        hull.tri[e] = a;
                        break;
                    }
                    e = hull.prev[e];
                    if e == hull.start {
                        break;
                    }
                }
            }

            self.halfedges[a] = hbl;
            self.halfedges[b] = har;
            self.halfedges[ar] = bl;

            if hbl != EMPTY {
                self.halfedges[hbl] = a;
            }
            if har != EMPTY {
                self.halfedges[har] = b;
            }
            if bl != EMPTY {
                self.halfedges[bl] = ar;
            }

            let br = next_halfedge(b);

            self.legalize(a, points, hull);
            return self.legalize(br, points, hull);
        }
        ar
    }
}

// data structure for tracking the edges of the advancing convex hull
struct Hull<T>
where
    T: CoordFloat,
{
    prev: Vec<usize>,
    next: Vec<usize>,
    tri: Vec<usize>,
    hash: Vec<usize>,
    start: usize,
    center: Point<T>,
}

impl<T> Hull<T>
where
    T: CoordFloat,
{
    fn new(
        n: usize,
        center: Point<T>,
        i0: usize,
        i1: usize,
        i2: usize,
        points: &[Point<T>],
    ) -> Self {
        let hash_len = (n as f64).sqrt() as usize;

        let mut hull = Self {
            prev: vec![0; n],            // edge to prev edge
            next: vec![0; n],            // edge to next edge
            tri: vec![0; n],             // edge to adjacent halfedge
            hash: vec![EMPTY; hash_len], // angular edge hash
            start: i0,
            center,
        };

        hull.next[i0] = i1;
        hull.prev[i2] = i1;
        hull.next[i1] = i2;
        hull.prev[i0] = i2;
        hull.next[i2] = i0;
        hull.prev[i1] = i0;

        hull.tri[i0] = 0;
        hull.tri[i1] = 1;
        hull.tri[i2] = 2;

        hull.hash_edge(points[i0], i0);
        hull.hash_edge(points[i1], i1);
        hull.hash_edge(points[i2], i2);

        hull
    }

    fn hash_key(&self, p: Point<T>) -> usize
    where
        T: CoordFloat,
    {
        let d = p - self.center;

        let p = d.x() / (d.x().abs() + d.y().abs());
        let a = (if d.y() > T::zero() {
            T::from(3.0).unwrap() - p
        } else {
            T::one() + p
        }) / T::from(4.0).unwrap(); // [0..1]

        let len = self.hash.len();
        let hash: T = NumCast::from(len).unwrap();
        let hash = (a * hash).floor();
        let hash: usize = NumCast::from(hash).unwrap();
        hash % len
    }

    fn hash_edge(&mut self, p: Point<T>, i: usize) {
        let key = self.hash_key(p);
        self.hash[key] = i;
    }

    fn find_visible_edge(&self, p: Point<T>, points: &[Point<T>]) -> (usize, bool) {
        let mut start: usize = 0;
        let key = self.hash_key(p);
        let len = self.hash.len();
        for j in 0..len {
            start = self.hash[(key + j) % len];
            if start != EMPTY && self.next[start] != EMPTY {
                break;
            }
        }
        start = self.prev[start];
        let mut e = start;

        while orient(p, points[e], points[self.next[e]]) {
            e = self.next[e];
            if e == start {
                return (EMPTY, false);
            }
        }
        (e, e == start)
    }
}

fn calc_bbox_center<T>(points: &[Point<T>]) -> Point<T>
where
    T: CoordFloat,
{
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

fn find_closest_point<T>(points: &[Point<T>], p0: Point<T>) -> Option<usize>
where
    T: CoordFloat,
{
    let mut min_dist = Float::infinity();
    let mut k: usize = 0;

    for (i, p) in points.iter().enumerate() {
        let d = dist2(p0, p);
        if d.is_sign_positive() && d < min_dist {
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

fn find_seed_triangle<T>(points: &[Point<T>]) -> Option<(usize, usize, usize)>
where
    T: CoordFloat,
{
    // pick a seed point close to the center
    let bbox_center = calc_bbox_center(points);
    let i0 = find_closest_point(points, bbox_center)?;
    let p0 = points[i0];

    // find the point closest to the seed
    let i1 = find_closest_point(points, p0)?;
    let p1 = points[i1];

    // find the third point which forms the smallest circumcircle with the first two
    let mut min_radius = Float::infinity();
    let mut i2: usize = 0;
    for (i, p) in points.iter().enumerate() {
        if i == i0 || i == i1 {
            continue;
        }
        let r = circumradius2(p0, p1, *p);
        if r < min_radius {
            i2 = i;
            min_radius = r;
        }
    }

    if min_radius == Float::infinity() {
        None
    } else {
        // swap the order of the seed points for counter-clockwise orientation
        Some(if orient(p0, p1, points[i2]) {
            (i0, i2, i1)
        } else {
            (i0, i1, i2)
        })
    }
}

fn sortf<T>(f: &mut Vec<(usize, T)>)
where
    T: CoordFloat,
{
    f.sort_unstable_by(|&(_, da), &(_, db)| da.partial_cmp(&db).unwrap());
}

/// Order collinear points by dx (or dy if all x are identical) and return the list as a hull
fn handle_collinear_points<T>(points: &[Point<T>]) -> Triangulation
where
    T: CoordFloat,
{
    let (x, y) = points
        .first()
        .map(|p| p.x_y())
        .unwrap_or((T::zero(), T::zero()));

    let mut dist: Vec<_> = points
        .iter()
        .enumerate()
        .map(|(i, p)| {
            let mut d = p.x() - x;
            if d.is_zero() {
                d = p.y() - y;
            }
            (i, d)
        })
        .collect();
    sortf(&mut dist);

    let mut triangulation = Triangulation::new(0);
    let mut d0 = Float::neg_infinity();
    for (i, distance) in dist {
        if distance > d0 {
            triangulation.hull.push(i);
            d0 = distance;
        }
    }

    triangulation
}

/// Triangulate a set of 2D points.
/// Returns the triangulation for the input points.
/// For the degenerated case when all points are collinear, returns an empty triangulation where all points are in the hull.
pub fn triangulate<T>(points: &[Point<T>]) -> Triangulation
where
    T: CoordFloat + AbsDiffEq<Epsilon = T>,
{
    let seed_triangle = find_seed_triangle(points);
    if seed_triangle.is_none() {
        return handle_collinear_points(points);
    }

    let n = points.len();
    let (i0, i1, i2) =
        seed_triangle.expect("At this stage, points are guaranteed to yeild a seed triangle");
    let center = circumcenter(points[i0], points[i1], points[i2]);

    let mut triangulation = Triangulation::new(n);
    triangulation.add_triangle(i0, i1, i2, EMPTY, EMPTY, EMPTY);

    // sort the points by distance from the seed triangle circumcenter
    let mut dists: Vec<_> = points
        .iter()
        .enumerate()
        .map(|(i, point)| (i, dist2(center, point)))
        .collect();

    sortf(&mut dists);

    let mut hull = Hull::new(n, center, i0, i1, i2, points);

    for (k, &(i, _)) in dists.iter().enumerate() {
        let p = points[i];

        // skip near-duplicates
        if k > 0 && abs_diff_eq!(p, points[dists[k - 1].0]) {
            continue;
        }
        // skip seed triangle points
        if i == i0 || i == i1 || i == i2 {
            continue;
        }

        // find a visible edge on the convex hull using edge hash
        let (mut e, walk_back) = hull.find_visible_edge(p, points);
        if e == EMPTY {
            continue; // likely a near-duplicate point; skip it
        }

        // add the first triangle from the point
        let t = triangulation.add_triangle(e, i, hull.next[e], EMPTY, EMPTY, hull.tri[e]);

        // recursively flip triangles from the point until they satisfy the Delaunay condition
        hull.tri[i] = triangulation.legalize(t + 2, points, &mut hull);
        hull.tri[e] = t; // keep track of boundary triangles on the hull

        // walk forward through the hull, adding more triangles and flipping recursively
        let mut n = hull.next[e];
        loop {
            let q = hull.next[n];
            if !orient(p, points[n], points[q]) {
                break;
            }
            let t = triangulation.add_triangle(n, i, q, hull.tri[i], EMPTY, hull.tri[n]);
            hull.tri[i] = triangulation.legalize(t + 2, points, &mut hull);
            hull.next[n] = EMPTY; // mark as removed
            n = q;
        }

        // walk backward from the other side, adding more triangles and flipping
        if walk_back {
            loop {
                let q = hull.prev[e];
                if !orient(p, points[q], points[e]) {
                    break;
                }
                let t = triangulation.add_triangle(q, i, e, EMPTY, hull.tri[e], hull.tri[q]);
                triangulation.legalize(t + 2, points, &mut hull);
                hull.tri[q] = t;
                hull.next[e] = EMPTY; // mark as removed
                e = q;
            }
        }

        // update the hull indices
        hull.prev[i] = e;
        hull.next[i] = n;
        hull.prev[n] = i;
        hull.next[e] = i;
        hull.start = e;

        // save the two new edges in the hash table
        hull.hash_edge(p, i);
        hull.hash_edge(points[e], e);
    }

    // expose hull as a vector of point indices
    let mut e = hull.start;
    loop {
        triangulation.hull.push(e);
        e = hull.next[e];
        if e == hull.start {
            break;
        }
    }

    triangulation.triangles.shrink_to_fit();
    triangulation.halfedges.shrink_to_fit();

    triangulation
}
