use crate::el_interpolation::calculate_witness_poly;
use crate::el_interpolation::calculate_zero_poly_coefficients;
use crate::el_interpolation::el_lagrange_interpolation;
use crate::el_interpolation::ElPoint;
use crate::toeplitz::ToeplitzMatrix;
use ark_bn254::{Bn254, Fr, G1Projective as G1, G2Projective as G2};
use ark_ec::pairing::Pairing;
use ark_ec::Group;
use ark_ff::Field;
use ark_poly::univariate::DensePolynomial;
use ark_poly::DenseUVPolynomial;
use ark_poly::Polynomial;
use ark_std::UniformRand;
use ark_std::Zero;
use rand::Rng;

// Entity that represents a random value that is calculated as a result of trusted setup.
// It could be generated using MPC
#[derive(Clone)]
pub struct CRS {
    pub powers_g1: Vec<G1>,
    pub powers_g2: Vec<G2>,
}

impl CRS {
    // FIXME - Insecure, should be used only for testing
    pub fn new(value: Fr, max_degree: usize) -> Self {
        let value_powers = calc_powers(value, max_degree);
        let powers_g1 = value_powers
            .iter()
            .map(|pow| G1::generator() * pow)
            .collect();
        let powers_g2 = value_powers
            .iter()
            .map(|pow| G2::generator() * pow)
            .collect();
        Self {
            powers_g1,
            powers_g2,
        }
    }

    // Returns a structure with random value
    // FIXME - Insecure, should be used only for testing
    pub fn new_rand<R: Rng>(rng: &mut R, max_degree: usize) -> Self {
        let value = Fr::rand(rng);
        let value_powers = calc_powers(value, max_degree);
        let powers_g1 = value_powers
            .iter()
            .map(|pow| G1::generator() * pow)
            .collect();
        let powers_g2 = value_powers
            .iter()
            .map(|pow| G2::generator() * pow)
            .collect();
        Self {
            powers_g1,
            powers_g2,
        }
    }

    // TODO: somewhere here an MPC protocol could be described or a Fiat-Shamir transform could be used
    pub fn mpc() {
        unimplemented!()
    }
}

// Calculates a powers of the value.
pub fn calc_powers(value: Fr, max_degree: usize) -> Vec<Fr> {
    let mut powers = Vec::with_capacity(max_degree + 1);
    for i in 0..=max_degree {
        powers.push(value.pow([i as u64]));
    }
    powers
}

#[derive(Debug)]
pub struct KZGProof {
    // I(X) - polynomial that passes through desired points for the check (zero at y)
    pub numerator: G1,
    // Z(X) - zero polynomial that has zeroes at xi (zeroes at x)
    pub denominator: G2,
    // q(x)
    pub witness: G1,
}
use crate::GeneralEvaluationDomain;
use ark_ff::FftField;
use ark_poly::EvaluationDomain;

impl KZGProof {
    pub fn new(numerator: G1, denominator: G2, witness: G1) -> Self {
        KZGProof {
            numerator,
            denominator,
            witness,
        }
    }

    // Prove function that follows the KZG procedure
    // FIXME: value: Fr is here only for testing purposes. Such a solution is not secure
    pub fn prove(
        _crs: &CRS,
        value: Fr,
        commit_poly: DensePolynomial<Fr>,
        witness_to: &[ElPoint],
    ) -> Self {
        // let powers = crs.calc_powers(commit_to.len());

        let i_coeffs = el_lagrange_interpolation(witness_to);
        let i_poly = DensePolynomial::from_coefficients_vec(i_coeffs);
        let numerator = G1::generator() * i_poly.evaluate(&value);

        // NOTE: I checked that this zero polynomial is generated correctly
        let zero_points: Vec<Fr> = witness_to.iter().map(|point| point.x).collect();
        let z_coeffs = calculate_zero_poly_coefficients(&zero_points);
        let z_poly = DensePolynomial::from_coefficients_vec(z_coeffs);
        let denominator = G2::generator() * z_poly.evaluate(&value);

        let witness_poly = calculate_witness_poly(&commit_poly, &i_poly, &z_poly);
        let witness = G1::generator() * witness_poly.evaluate(&value);

        KZGProof {
            numerator,
            denominator,
            witness,
        }
    }

    // The proof verification for one point: e(q(x)1, [x-x']2) == e([p(x)-p(x')]1, G2)
    // The proof verification for several points: e(q(x)1, [Z(x)]2) == e([p(x)-I(x)]1, G2)
    //FIXME - Insecure, there should be check that the numerator and denumerator are calculated correctly
    pub fn verify(&self, commitment: G1) -> bool {
        let left = Bn254::pairing(self.witness, self.denominator);
        let right = Bn254::pairing(commitment - self.numerator, G2::generator());

        left == right
    }

    pub fn calc_all_proofpoints(crs: &CRS, commit_to: &[ElPoint]) -> Vec<G1> {
        let mut commit_coeffs = el_lagrange_interpolation(commit_to);
        commit_coeffs.reverse();

        // a vector that has the first coefficient from the commit_coeffs and all other values are Fr::zero()
        let mut zeros = vec![Fr::zero(); commit_coeffs.len()];
        zeros[0] = commit_coeffs[0];
        let toeplitz = ToeplitzMatrix::new(commit_coeffs, zeros).unwrap();
        let circulant = toeplitz.extend_to_circulant();
        let hs = circulant.fast_multiply_by_vector(&crs.powers_g1).unwrap();
        hs
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_std::One;
    use ark_std::Zero;

    #[test]
    fn test_fr_calc_powers() {
        let value = Fr::from(1423);
        let max_degree: usize = 4;
        let crs = CRS::new(value, max_degree);

        assert_eq!(crs.powers_g1.len(), max_degree + 1);
        assert_eq!(crs.powers_g2.len(), max_degree + 1);

        assert_eq!(crs.powers_g1[0], G1::generator() * Fr::one());
        assert_eq!(crs.powers_g1[1], G1::generator() * value);
        assert_eq!(crs.powers_g1[2], G1::generator() * value * value);
        assert_eq!(crs.powers_g1[3], G1::generator() * value * value * value);
        assert_eq!(
            crs.powers_g1[4],
            G1::generator() * value * value * value * value
        );
    }

    #[test]
    fn test_kzg_proof_verify_valid() {
        let commit_to = vec![
            ElPoint::new(Fr::from(1), Fr::from(2)),
            ElPoint::new(Fr::from(2), Fr::from(3)),
            ElPoint::new(Fr::from(3), Fr::from(5)),
        ];
        let witness_to = vec![
            ElPoint::new(Fr::from(1), Fr::from(2)),
            ElPoint::new(Fr::from(2), Fr::from(3)),
        ];

        let max_degree = commit_to.len();
        let value = Fr::from(1423);
        let crs = CRS::new(value, max_degree);

        let commit_coeffs = el_lagrange_interpolation(&commit_to);
        let commit_poly = DensePolynomial::from_coefficients_vec(commit_coeffs);
        let commitment = G1::generator() * commit_poly.evaluate(&value);

        let proof = KZGProof::prove(&crs, value, commit_poly, &witness_to);

        assert!(proof.verify(commitment));
    }

    #[test]
    fn test_kzg_proof_verify_valid_commit_and_witness_equals() {
        let commit_to = vec![
            ElPoint::new(Fr::from(1), Fr::from(2)),
            ElPoint::new(Fr::from(2), Fr::from(3)),
            ElPoint::new(Fr::from(3), Fr::from(5)),
        ];
        let witness_to = vec![
            ElPoint::new(Fr::from(1), Fr::from(2)),
            ElPoint::new(Fr::from(2), Fr::from(3)),
            ElPoint::new(Fr::from(3), Fr::from(5)),
        ];

        let max_degree = commit_to.len();
        let value = Fr::from(1423);
        let crs = CRS::new(value, max_degree);

        let commit_coeffs = el_lagrange_interpolation(&commit_to);
        let commit_poly = DensePolynomial::from_coefficients_vec(commit_coeffs);
        let commitment = G1::generator() * commit_poly.evaluate(&value);

        let proof = KZGProof::prove(&crs, value, commit_poly, &witness_to);

        assert!(proof.verify(commitment));
    }

    #[test]
    fn test_kzg_proof_verify_valid_zero() {
        let commitment = G1::zero();
        let numerator = G1::zero();
        let denominator = G2::zero();
        let witness = G1::zero();

        let proof = KZGProof {
            numerator,
            denominator,
            witness,
        };

        assert!(proof.verify(commitment));
    }

    #[test]
    fn test_kzg_proof_verify_invalid() {
        let mut rng: rand::prelude::ThreadRng = rand::thread_rng();
        let commitment = G1::rand(&mut rng);
        let numerator = G1::rand(&mut rng);
        let denominator = G2::rand(&mut rng);
        let witness = G1::rand(&mut rng);

        let proof = KZGProof {
            numerator,
            denominator,
            witness,
        };

        assert!(!proof.verify(commitment));
    }

    #[test]
    fn test_calc_all_proofpoints() {
        // Sample CRS
        let crs = CRS::new(Fr::from(6), 5);

        // Sample commit_to vector
        let commit_to = vec![
            ElPoint::new(Fr::from(1), Fr::from(2)),
            ElPoint::new(Fr::from(2), Fr::from(3)),
            ElPoint::new(Fr::from(3), Fr::from(5)),
        ];

        // Calculate the proof points
        let proof_points = KZGProof::calc_all_proofpoints(&crs, &commit_to);

        println!("{:#?}", proof_points);
    }
}
