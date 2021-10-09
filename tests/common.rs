use std::fs::File;

use delaunator::{triangulate, CoordType, EMPTY};
use geo_types::{point, CoordFloat, Point};

fn validate<T>(points: &[Point<T>])
where
    T: CoordType,
{
    let triangulation = triangulate(&points);

    // validate halfedges
    for (i, &h) in triangulation.halfedges.iter().enumerate() {
        if h != EMPTY && triangulation.halfedges[h] != i {
            panic!("Invalid halfedge connection");
        }
    }

    // validate triangulation
    let hull_area = {
        let mut hull_areas = Vec::new();
        let mut i = 0;
        let mut j = triangulation.hull.len() - 1;
        while i < triangulation.hull.len() {
            let p0 = &points[triangulation.hull[j]];
            let p = &points[triangulation.hull[i]];
            hull_areas.push((p.x() - p0.x()) * (p.y() + p0.y()));
            j = i;
            i += 1;
        }
        sum(&hull_areas)
    };
    let triangles_area = {
        let mut triangle_areas = Vec::new();
        let mut i = 0;
        while i < triangulation.triangles.len() {
            let a = &points[triangulation.triangles[i]];
            let b = &points[triangulation.triangles[i + 1]];
            let c = &points[triangulation.triangles[i + 2]];
            triangle_areas.push(
                ((b.y() - a.y()) * (c.x() - b.x()) - (b.x() - a.x()) * (c.y() - b.y())).abs(),
            );
            i += 3;
        }
        sum(&triangle_areas)
    };

    let err = ((hull_area - triangles_area) / hull_area).abs();
    if err > T::epsilon() {
        panic!("Triangulation is broken: {} error", err);
    }
}

pub fn basic<T>()
where
    T: CoordType,
{
    validate::<T>(&load_fixture("tests/fixtures/basic.json"));
}

pub fn js_issues<T>()
where
    T: CoordType,
{
    validate::<T>(&load_fixture("tests/fixtures/issue11.json"));
    validate::<T>(&load_fixture("tests/fixtures/issue13.json"));
    validate::<T>(&load_fixture("tests/fixtures/issue24.json"));
}

pub fn robustness<T>()
where
    T: CoordType + CoordFloat,
{
    let points = load_fixture::<T>("tests/fixtures/robust1.json");

    validate::<T>(&points);
    validate::<T>(&(scale_points(&points, T::from(1e-9).unwrap())));
    validate::<T>(&(scale_points(&points, T::from(1e-2).unwrap())));
    validate::<T>(&(scale_points(&points, T::from(100.0).unwrap())));
    validate::<T>(&(scale_points(&points, T::from(1e9).unwrap())));

    validate::<T>(&load_fixture("tests/fixtures/robust2.json"));
    validate::<T>(&load_fixture("tests/fixtures/robust3.json"));
}

pub fn bad_input<T>()
where
    T: CoordType,
{
    let mut points = vec![];
    let triangulation = triangulate(&points);

    assert!(
        triangulation.triangles.is_empty(),
        "Expected no triangles (0 point)"
    );
    assert!(
        triangulation.halfedges.is_empty(),
        "Expected no edges (0 point)"
    );
    assert!(triangulation.hull.is_empty(), "Expected no hull (0 point)");

    points.push(point!(x: 0., y: 0. ));
    let triangulation = triangulate(&points);

    assert!(
        triangulation.triangles.is_empty(),
        "Expected no triangles (1 point)"
    );
    assert!(
        triangulation.halfedges.is_empty(),
        "Expected no edges (1 point)"
    );
    assert!(
        triangulation.hull.len() == 1,
        "Expected single point on hull (1 point)"
    );

    points.push(point!(x: 1., y: 0. ));
    let triangulation = triangulate(&points);

    assert!(
        triangulation.triangles.is_empty(),
        "Expected no triangles (2 points)"
    );
    assert!(
        triangulation.halfedges.is_empty(),
        "Expected no edges (2 points)"
    );
    assert!(
        triangulation.hull.len() == 2,
        "Expected two points on hull (2 point)"
    );
    assert!(
        triangulation.hull.iter().enumerate().all(|(i, v)| i == *v),
        "Expected ordered hull points (2 point)"
    );

    points.push(point!(x: 2., y: 0. ));
    let triangulation = triangulate(&points);

    assert!(
        triangulation.triangles.is_empty(),
        "Expected no triangles (3 collinear points)"
    );
    assert!(
        triangulation.halfedges.is_empty(),
        "Expected no edges (3 collinear points)"
    );
    assert!(
        triangulation.hull.len() == 3,
        "Expected three points on hull (3 collinear points)"
    );
    assert!(
        triangulation.hull.iter().enumerate().all(|(i, v)| i == *v),
        "Expected ordered hull points (3 collinear points)"
    );

    points.push(point!(x: 1., y: 1. ));
    validate(&points);
}

pub fn unordered_collinear_points_input<T>()
where
    T: CoordType,
{
    let points: Vec<Point<T>> = [10, 2, 4, 4, 1, 0, 3, 6, 8, 5, 7, 9]
        .iter()
        .map(|i| T::from(*i).unwrap())
        .map(|y| point!(x: T::zero(), y: y))
        .collect();
    let duplicated = 1;

    let triangulation = triangulate(&points);

    assert!(
        triangulation.triangles.is_empty(),
        "Expected no triangles (unordered collinear points)"
    );
    assert!(
        triangulation.halfedges.is_empty(),
        "Expected no edges (unordered collinear points)"
    );
    assert!(
        triangulation.hull.len() == points.len() - duplicated,
        "Expected all non-coincident points on hull (unordered collinear points)"
    );
    assert!(
        triangulation
            .hull
            .iter()
            .enumerate()
            .all(|(i, v)| points[*v].y() == T::from(i).unwrap()),
        "Expected ordered hull points (unordered collinear points)"
    );
}

fn scale_points<T>(points: &[Point<T>], scale: T) -> Vec<Point<T>>
where
    T: CoordType,
{
    let scaled = points
        .iter()
        .map(|p| {
            point!(
                x: p.x() * scale,
                y: p.y() * scale
            )
        })
        .collect();
    scaled
}

fn load_fixture<T>(path: &str) -> Vec<Point<T>>
where
    T: CoordType,
{
    let file = File::open(path).unwrap();
    let u: Vec<(f64, f64)> = serde_json::from_reader(file).unwrap();
    u.iter()
        .map(|p| {
            let x = T::from(p.0);
            assert!(
                x.is_some(),
                "Unable to parse {} into {}",
                p.0,
                std::any::type_name::<T>()
            );

            let y = T::from(p.1);
            assert!(
                y.is_some(),
                "Unable to parse {} into {}",
                p.1,
                std::any::type_name::<T>()
            );
            point!(x: x.unwrap(), y: y.unwrap())
        })
        .collect()
}

// Kahan and Babuska summation, Neumaier variant; accumulates less FP error
fn sum<T>(x: &[T]) -> T
where
    T: CoordType,
{
    let mut sum = x[0];
    let mut err: T = T::zero();
    for i in 1..x.len() {
        let k = x[i];
        let m = sum + k;
        err = err
            + if sum.abs() >= k.abs() {
                sum - m + k
            } else {
                k - m + sum
            };
        sum = m;
    }
    sum + err
}
