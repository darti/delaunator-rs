use delaunator::{triangulate, Triangulation};
use geo_types::{point, Point};
use iai::{black_box, main};
use rand::{rngs::StdRng, Rng, SeedableRng};
use std::iter::repeat_with;

fn bench(count: usize) -> Triangulation<f64> {
    let mut rng: StdRng = StdRng::seed_from_u64(123);

    let points: Vec<Point<f64>> = repeat_with(|| rng.gen())
        .map(|(x, y)| point!(x: x, y: y))
        .take(count)
        .collect();

    triangulate(&points)
}

fn iai_bench_100() -> Triangulation<f64> {
    bench(black_box(100))
}

fn iai_bench_1000() -> Triangulation<f64> {
    bench(black_box(1_000))
}

fn iai_bench_10000() -> Triangulation<f64> {
    bench(black_box(10_000))
}

fn iai_bench_100000() -> Triangulation<f64> {
    bench(black_box(100_000))
}

iai::main!(
    iai_bench_100,
    iai_bench_1000,
    iai_bench_10000,
    iai_bench_100000
);
