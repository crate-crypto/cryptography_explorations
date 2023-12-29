pub mod el_interpolation;
pub mod group;
pub mod kzg;
pub mod kzg_transcript;
pub mod pedersen;
pub mod toeplitz;

use crate::el_interpolation::el_lagrange_interpolation;
use crate::el_interpolation::ElPoint;
use ark_bn254::Fr;
use ark_poly::univariate::DensePolynomial;
use ark_poly::DenseUVPolynomial;
use ark_poly::EvaluationDomain;
use ark_poly::GeneralEvaluationDomain;
use std::ops::Mul;

// use nalgebra::base::Matrix4;
// use nalgebra::base::Vector4;
fn main() {
    // Experiments with FFT
    let poly_evals1 = vec![
        ElPoint::new(Fr::from(1), Fr::from(2)),
        ElPoint::new(Fr::from(2), Fr::from(3)),
        ElPoint::new(Fr::from(3), Fr::from(5)),
    ];
    let poly_evals2 = vec![
        ElPoint::new(Fr::from(1), Fr::from(7)),
        ElPoint::new(Fr::from(2), Fr::from(9)),
    ];
    let coeff1 = el_lagrange_interpolation(&poly_evals1);
    let coeff2 = el_lagrange_interpolation(&poly_evals2);

    let domain: GeneralEvaluationDomain<Fr> = GeneralEvaluationDomain::new(4).unwrap();
    let fft1 = domain.fft(&coeff1);
    let fft2 = domain.fft(&coeff2);

    let mul_poly_evals = domain.mul_polynomials_in_evaluation_domain(&fft1, &fft2);
    let mul_coeffs = domain.ifft(&mul_poly_evals);
    println!("{:#?}", mul_coeffs);

    // O(n log n) multiplication
    let poly1 = DensePolynomial::from_coefficients_vec(coeff1.to_vec());
    let poly2 = DensePolynomial::from_coefficients_vec(coeff2.to_vec());
    let mul1 = poly1.mul(&poly2);
    print!("{:#?}", mul1);

    // let vec2 = Vector4::new(Fr::from(5), Fr::from(6), Fr::from(7), Fr::from(8));
    // // Multiplication happens
    // let mat = Matrix4::new(
    //     Fr::from(1),
    //     Fr::from(4),
    //     Fr::from(3),
    //     Fr::from(2),
    //     Fr::from(2),
    //     Fr::from(1),
    //     Fr::from(4),
    //     Fr::from(3),
    //     Fr::from(3),
    //     Fr::from(2),
    //     Fr::from(1),
    //     Fr::from(4),
    //     Fr::from(4),
    //     Fr::from(3),
    //     Fr::from(2),
    //     Fr::from(1),
    // );

    // let mul = mat * vec2;
    // println!(
    //     "{:?}",
    //     mul.iter().map(|elem| elem.to_string()).collect::<Vec<_>>()
    // );
}
