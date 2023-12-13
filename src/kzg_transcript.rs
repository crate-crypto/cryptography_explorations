use crate::el_interpolation::calculate_witness_poly;
use crate::el_interpolation::calculate_zero_poly_coefficients;
use crate::el_interpolation::el_lagrange_interpolation;
use crate::el_interpolation::ElPoint;
use ark_bn254::{Bn254, Fr, G1Projective as G1, G2Projective as G2};
use ark_ec::pairing::Pairing;
use ark_ec::Group;
use ark_ff::Field;
use ark_poly::univariate::DensePolynomial;
use ark_poly::DenseUVPolynomial;
use ark_poly::Polynomial;
use ark_std::UniformRand;
use rand::Rng;

// TODO: GPU Fast Fourier Transform and Multiexponentiation:
// https://docs.rs/ec-gpu-gen/latest/ec_gpu_gen/
// https://github.com/filecoin-project/bellperson/blob/master/src/gpu/multiexp.rs

// Entity that represents a random value that is calculated as a result of trusted setup.
// It could be generated using MPC
//FIXME - Insecure, should be used only for testing
#[derive(Copy, Clone)]
pub struct CRS {
    pub value: Fr,
}

impl CRS {
    // TODO: somewhere here an MPC protocol could be described or a Fiat-Shamir transform could be used
    //FIXME - Insecure, should be used only for testing
    pub fn new(value: Fr) -> Self {
        Self { value }
    }

    // Returns a structure with random value
    //FIXME - Insecure, should be used only for testing
    pub fn new_rand<R: Rng>(rng: &mut R) -> Self {
        Self {
            value: Fr::rand(rng),
        }
    }

    // Calculates a powers of the value.
    //FIXME - Insecure, should be used only for testing
    pub fn calc_powers(&self, max_degree: usize) -> Vec<Fr> {
        let mut powers = Vec::with_capacity(max_degree + 1);
        for i in 0..=max_degree {
            powers.push(self.value.pow([i as u64]));
        }
        powers
    }
}

pub struct KZGParams {
    pub max_degree: usize,
    pub powers: Vec<Fr>,
}

impl KZGParams {
    pub fn new(max_degree: usize, powers: Vec<Fr>) -> Self {
        Self { max_degree, powers }
    }
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

impl KZGProof {
    pub fn new(numerator: G1, denominator: G2, witness: G1) -> Self {
        KZGProof {
            numerator,
            denominator,
            witness,
        }
    }

    // Prove function that follows the KZG procedure
    pub fn prove(crs: &CRS, commit_poly: DensePolynomial<Fr>, witness_to: &[ElPoint]) -> Self {
        // let powers = crs.calc_powers(commit_to.len());

        let i_coeffs = el_lagrange_interpolation(witness_to);
        let i_poly = DensePolynomial::from_coefficients_vec(i_coeffs);
        let numerator = G1::generator() * i_poly.evaluate(&crs.value);

        // NOTE: I checked that this zero polynomial is generated correctly
        let zero_points: Vec<Fr> = witness_to.iter().map(|point| point.x).collect();
        let z_coeffs = calculate_zero_poly_coefficients(&zero_points);
        let z_poly = DensePolynomial::from_coefficients_vec(z_coeffs);
        let denominator = G2::generator() * z_poly.evaluate(&crs.value);

        let witness_poly = calculate_witness_poly(&commit_poly, &i_poly, &z_poly);
        let witness = G1::generator() * witness_poly.evaluate(&crs.value);

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
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_std::One;
    use ark_std::Zero;

    #[test]
    fn test_fr_calc_powers() {
        let value: Fr = Fr::from(2);
        let crs = CRS::new(value);
        let max_degree: usize = 4;
        let powers = crs.calc_powers(max_degree);

        assert_eq!(powers.len(), max_degree + 1);
        assert_eq!(powers[0], Fr::one());
        assert_eq!(powers[1], value);
        assert_eq!(powers[2], value * value);
        assert_eq!(powers[3], value * value * value);
        assert_eq!(powers[4], value * value * value * value);
    }

    #[test]
    fn test_kzg_proof_verify_valid() {
        let crs = CRS::new(Fr::from(5));

        let commit_to = vec![
            ElPoint::new(Fr::from(1), Fr::from(2)),
            ElPoint::new(Fr::from(2), Fr::from(3)),
            ElPoint::new(Fr::from(3), Fr::from(5)),
        ];
        let witness_to = vec![
            ElPoint::new(Fr::from(1), Fr::from(2)),
            ElPoint::new(Fr::from(2), Fr::from(3)),
        ];

        let commit_coeffs = el_lagrange_interpolation(&commit_to);
        let commit_poly = DensePolynomial::from_coefficients_vec(commit_coeffs);
        let commitment = G1::generator() * commit_poly.evaluate(&crs.value);

        let proof = KZGProof::prove(&crs, commit_poly, &witness_to);

        assert!(proof.verify(commitment));
    }

    #[test]
    fn test_kzg_proof_verify_valid_commit_and_witness_equals() {
        let crs = CRS::new(Fr::from(13));

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

        let commit_coeffs = el_lagrange_interpolation(&commit_to);
        let commit_poly = DensePolynomial::from_coefficients_vec(commit_coeffs);
        let commitment = G1::generator() * commit_poly.evaluate(&crs.value);

        let proof = KZGProof::prove(&crs, commit_poly, &witness_to);

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
}
