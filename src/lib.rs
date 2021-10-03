use approx::AbsDiffEq;
pub use delaunay::{Triangulation, EMPTY};
use geo_types::{CoordFloat, Point};
use math::CoordType;

mod delaunay;
mod math;

pub fn triangulate<T>(points: &[Point<T>]) -> Triangulation<T>
where
    T: CoordType,
{
    Triangulation::triangulate(points)
}
