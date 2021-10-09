pub use delaunay::Triangulation;
use geo_types::Point;
pub use math::{CoordType, EMPTY};

mod delaunay;
mod math;

pub fn triangulate<T>(points: &[Point<T>]) -> Triangulation<T>
where
    T: CoordType,
{
    Triangulation::triangulate(points)
}
