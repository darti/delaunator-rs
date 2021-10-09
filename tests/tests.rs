use common::{bad_input, basic, js_issues, robustness, unordered_collinear_points_input};
use paste::paste;

mod common;

macro_rules! test_type {
    ($name: ident, $($typ: ty), +) => {
        $(
        paste! {
            #[test]
            fn [<$name _ $typ>]() {
                $name::<$typ>();
            }
        })+
    };
}

test_type!(basic, f64, f32, i64);
test_type!(js_issues, f64, f32);
test_type!(robustness, f64);
test_type!(bad_input, f64, f32);
test_type!(unordered_collinear_points_input, f64, f32);
