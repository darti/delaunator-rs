use approx::AbsDiffEq;
pub use delaunay::Triangulation;
use geo_types::{CoordFloat, Point};

mod delaunay;

pub fn triangulate<T>(points: &[Point<T>]) -> Triangulation<T>
where
    T: CoordFloat + AbsDiffEq<Epsilon = T>,
{
    Triangulation::triangulate(points)
}
