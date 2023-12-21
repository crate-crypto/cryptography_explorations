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

use nalgebra::base::Matrix4;
use nalgebra::base::Vector4;
fn main() {
    // Experiments with FFT
    let commit_to = vec![
        ElPoint::new(Fr::from(1), Fr::from(2)),
        ElPoint::new(Fr::from(2), Fr::from(3)),
        ElPoint::new(Fr::from(3), Fr::from(5)),
    ];
    let commit_to1 = vec![
        ElPoint::new(Fr::from(1), Fr::from(7)),
        ElPoint::new(Fr::from(2), Fr::from(9)),
    ];
    let commit_coeff = el_lagrange_interpolation(&commit_to);
    let commit_coeff1 = el_lagrange_interpolation(&commit_to1);

    let domain: GeneralEvaluationDomain<Fr> = GeneralEvaluationDomain::new(6).unwrap();
    let fft = domain.fft(&commit_coeff);
    let fft1 = domain.fft(&commit_coeff1);
    // println!("FFT {:#?}", fft);

    let mul = domain.mul_polynomials_in_evaluation_domain(&fft, &fft1);
    let mul_coeffs = domain.ifft(&mul);
    // println!("{:#?}", mul_coeffs);

    let commit_poly = DensePolynomial::from_coefficients_vec(commit_coeff.to_vec());
    let commit_poly1 = DensePolynomial::from_coefficients_vec(commit_coeff1.to_vec());
    let mul1 = commit_poly.mul(&commit_poly1);

    // print!("{:#?}", mul1);

    let vec2 = Vector4::new(Fr::from(5), Fr::from(6), Fr::from(7), Fr::from(8));
    // Multiplication happens 
    let mat = Matrix4::new(
        Fr::from(1),
        Fr::from(4),
        Fr::from(3),
        Fr::from(2),
        Fr::from(2),
        Fr::from(1),
        Fr::from(4),
        Fr::from(3),
        Fr::from(3),
        Fr::from(2),
        Fr::from(1),
        Fr::from(4),
        Fr::from(4),
        Fr::from(3),
        Fr::from(2),
        Fr::from(1),
    );

    // got: Ok([BigInt([66, 0, 0, 0]), BigInt([68, 0, 0, 0]), BigInt([66, 0, 0, 0]), BigInt([60, 0, 0, 0])])

    // Fr::from(1), Fr::from(2), Fr::from(3), Fr::from(4),
    //  Fr::from(2), Fr::from(3), Fr::from(4), Fr::from(1),
    //  Fr::from(3), Fr::from(4),Fr::from(1), Fr::from(2),
    //  Fr::from(4),Fr::from(1), Fr::from(2), Fr::from(3),

    let mul = mat * vec2;
    println!(
        "{:?}",
        mul.iter().map(|elem| elem.to_string()).collect::<Vec<_>>()
    );
}
