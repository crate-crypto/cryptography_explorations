// use crate::el_interpolation::calculate_zero_poly_coefficients;
// use crate::el_interpolation::{el_lagrange_interpolation, ElPoint};
// use crate::kzg_transcript::KZGParams;
// use crate::kzg_transcript::KZGProof;
// use crate::kzg_transcript::CRS;
// use ark_bn254::{Fr, G1Affine as G1, G2Affine as G2};
// use ark_ec::{AffineRepr, CurveGroup};
// use ark_poly::univariate::DensePolynomial;
// use ark_poly::DenseUVPolynomial;
// use ark_poly::Polynomial;

// #[derive(Debug, Clone, Copy)]
// pub struct KZGCommitment {
//     pub point: G1,
// }

// impl KZGCommitment {
//     pub fn setup(value: CRS, max_degree: usize) -> KZGParams {
//         KZGParams::new(max_degree, value.calc_powers(max_degree))
//     }

//     // e(q(x)1, [Z(x)]2) == e([p(x)-I(x)]1, G2)
//     // numerator - I(X),
//     // denominator - Z(X),
//     // witness - q(x),
//     // Proof: q(x) = (p(x) - I(x)) / Z(x)
//     pub fn commit_poly(coefficients: &[Fr], point: CRS, points: &[ElPoint]) -> (Self, KZGProof) {
//         let polynomial = DensePolynomial::from_coefficients_vec(coefficients.to_vec());
//         let commitment_raw = polynomial.evaluate(&point.value);
//         let commitment: ark_ec::short_weierstrass::Affine<ark_bn254::g1::Config> =
//             (G1::generator() * commitment_raw).into();
//         println!("commitment: {:?}", commitment_raw.to_string());

//         let numerator_coeffs = el_lagrange_interpolation(points);
//         let numerator_poly = DensePolynomial::from_coefficients_vec(numerator_coeffs.to_vec());
//         let numerator_raw = numerator_poly.evaluate(&point.value);
//         let numerator: ark_ec::short_weierstrass::Affine<ark_bn254::g1::Config> =
//             (G1::generator() * numerator_raw).into();
//         println!("numerator: {:?}", numerator_raw.to_string());

//         let zero_points: Vec<Fr> = points.iter().map(|point| point.x).collect();
//         let denominator_coeffs = calculate_zero_poly_coefficients(&zero_points);
//         let denominator_poly = DensePolynomial::from_coefficients_vec(denominator_coeffs.to_vec());
//         let denominator_raw = denominator_poly.evaluate(&point.value);
//         let denominator: ark_ec::short_weierstrass::Affine<ark_bn254::g2::Config> =
//             (G2::generator() * denominator_raw).into();
//         println!("denominator: {:?}", denominator_raw.to_string());

//         let witness_raw = (commitment_raw - numerator_raw) / denominator_raw;
//         let witness: ark_ec::short_weierstrass::Affine<ark_bn254::g1::Config> =
//             (G1::generator() * witness_raw).into();
//         let proof = KZGProof::prove(numerator, denominator, witness);
//         println!("witness: {:?}", witness_raw.to_string());

//         (KZGCommitment { point: commitment }, proof)
//     }

//     pub fn verify_poly(commitment: Self, coefficients: &[G1], value_powers: &Vec<Fr>) -> bool {
//         let mut res: ark_ec::short_weierstrass::Affine<ark_bn254::g1::Config> = G1::zero();

//         for (i, coefficient) in coefficients.iter().enumerate() {
//             let power = value_powers[i];
//             let result: ark_ec::short_weierstrass::Affine<ark_bn254::g1::Config> =
//                 (*coefficient * &power).into_affine();
//             res = (res + result).into_affine()
//         }

//         commitment.point == res
//     }

//     pub fn verify_witness(commitment: Self, proof: KZGProof) -> bool {
//         proof.verify(commitment.point)
//     }
// }

// #[cfg(test)]
// mod tests {
//     use super::*;
//     use rand::Rng;

//     #[test]
//     fn full_cycle_test() {
//         let point = CRS::new(Fr::from(13));

//         let points = vec![
//             ElPoint::new(Fr::from(1), Fr::from(2)),
//             ElPoint::new(Fr::from(2), Fr::from(3)),
//             ElPoint::new(Fr::from(3), Fr::from(5)),
//         ];
//         let points_for_check = vec![ElPoint::new(Fr::from(1), Fr::from(9))];

//         let commitment_coeff = el_lagrange_interpolation(&points);

//         let (commitment, proof) =
//             KZGCommitment::commit_poly(&commitment_coeff, point.clone(), &points_for_check);

//         // Verifier part
//         let numerator_coeffs = el_lagrange_interpolation(&points_for_check);
//         let numerator_poly = DensePolynomial::from_coefficients_vec(numerator_coeffs.to_vec());
//         let numerator_raw = numerator_poly.evaluate(&point.value);
//         println!("numerator_raw {:?}", numerator_raw.to_string());
//         let numerator: ark_ec::short_weierstrass::Affine<ark_bn254::g1::Config> =
//             (G1::generator() * numerator_raw).into();

//         let zero_points: Vec<Fr> = points_for_check.iter().map(|point| point.x).collect();
//         let denominator_coeffs = calculate_zero_poly_coefficients(&zero_points);
//         let denominator_poly = DensePolynomial::from_coefficients_vec(denominator_coeffs.to_vec());
//         let denominator_raw = denominator_poly.evaluate(&point.value);
//         println!("denominator_raw {:?}", denominator_raw.to_string());
//         let denominator: ark_ec::short_weierstrass::Affine<ark_bn254::g2::Config> =
//             (G2::generator() * denominator_raw).into();

//         let verifier_proof = KZGProof::new(commitment, numerator, denominator, proof.witness);
//         assert!(KZGCommitment::verify_witness(commitment, verifier_proof));
//     }

//     #[test]
//     fn test_setup() {
//         let mut rng: rand::prelude::ThreadRng = rand::thread_rng();
//         let value = CRS::new_rand(&mut rng);
//         // Generate a rand max_degree
//         let max_degree = rng.gen_range(1..10);
//         let params = KZGCommitment::setup(value, max_degree);
//         assert!(params.max_degree == max_degree);
//         assert!(params.powers.len() == max_degree + 1);
//     }
// }
