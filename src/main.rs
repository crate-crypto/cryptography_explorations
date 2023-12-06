use ark_bn254::G1Projective as G1;

use ark_poly::univariate::DensePolynomial;
use ark_poly::DenseUVPolynomial;

use ark_poly::Polynomial;
use ark_std::UniformRand;

pub mod el_interpolation;
use crate::el_interpolation::{el_lagrange_interpolation, ElPoint};

fn main() {
    let mut rng = rand::thread_rng();

    let num_points = 3;
    let mut el_points = Vec::new();

    // Generate the defined quantity of points
    for _ in 0..num_points {
        let point = G1::rand(&mut rng);
        el_points.push(ElPoint::new(point.x, point.y));
    }
    println!("Elliptic curve points: {el_points:#?}");

    // Calculate the polynomial coefficients using interpolation
    let coefficients = el_lagrange_interpolation(&el_points);
    let str_coeff: Vec<String> = coefficients.iter().map(|p| p.to_string()).collect();
    println!("Coefficients: {str_coeff:#?}");

    // Cunstruct a polynomial with defined coefficients
    let el_poly = DensePolynomial::from_coefficients_vec(coefficients);
    println!("Polynomial that passes through elliptic curve points: {el_poly:#?}");

    for point in el_points {
        println!(
            "Evaluation at point X: {:#?}: \n Y: {:#?}",
            point.x,
            el_poly.evaluate(&point.x)
        );
    }
}
