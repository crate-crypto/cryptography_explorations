pub mod el_interpolation;
pub mod group;
pub mod interpolation;
pub mod kzg;
pub mod kzg_transcript;
pub mod pedersen;
use ark_bn254::g1::G1Affine;
use ark_bn254::Fr;
use ark_ec::AffineRepr;
use ark_poly::univariate::DensePolynomial;
use ark_poly::univariate::SparsePolynomial;
use ark_test_curves::fp128::Fq;
use std::ops::SubAssign;

fn main() {
    // 12160134928799597345692447636254041715860202444675574635387891214764338003182 =
    // = (47 - 244444) / 9
    let a = Fr::from(47);
    let b = Fr::from(449494);
    let c = Fr::from(9);
    println!("{:?}", ((a - b) / c).to_string());

    
}
