pub mod el_interpolation;
pub mod group;
pub mod interpolation;
pub mod kzg;
pub mod kzg_transcript;
pub mod pedersen;
use crate::el_interpolation::*;
use ark_bn254::Fr;

fn main() {
    // Example usage:
    let points = vec![
        ElPoint {
            x: Fr::from(0),
            y: Fr::from(0),
        },
        ElPoint {
            x: Fr::from(2),
            y: Fr::from(4),
        },
        ElPoint {
            x: Fr::from(4),
            y: Fr::from(16),
        },
    ];

    let x_coordinate = Fr::from(6);

    //
    let result_y = barycentric_interpolation(&points, x_coordinate);

    println!("Interpolated Y coordinate: {:?}", result_y);
}
