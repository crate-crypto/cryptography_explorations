use ark_bn254::Bn254;
use ark_bn254::{Fr, G1Affine as G1, G2Affine as G2};
use ark_ec::pairing::Pairing;
use ark_ec::AffineRepr;
use ark_ff::Field;
use ark_std::UniformRand;
use rand::Rng;

// Entity that represents a random value that is calculated as a result of trusted setup.
// It chould be generated using MPC
#[derive(Copy, Clone)]
pub struct CRS {
    pub value: Fr,
}

impl CRS {
    // TODO: somewhere here an MPC protocol could be described or a Fiat-Shamir transform could be used
    pub fn new(value: Fr) -> Self {
        Self { value }
    }

    // Returns a structure with random value
    pub fn new_rand<R: Rng>(rng: &mut R) -> Self {
        Self {
            value: Fr::rand(rng),
        }
    }

    // Calculates a powers of the value
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

    // TODO
    pub fn prove(numerator: G1, denominator: G2, witness: G1) -> Self {
        KZGProof {
            numerator,
            denominator,
            witness,
        }
    }

    // The proof verification for one point: e(q(x)1, [x-x']2) == e([p(x)-p(x')]1, G2)
    // The proof verification for several points: e(q(x)1, [Z(x)]2) == e([p(x)-I(x)]1, G2)
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
    // TODO: make up a test case
    fn test_kzg_proof_verify_valid() {
        // let commitment = G1::rand(&mut rng);
        // let numerator = G1::rand(&mut rng);
        // let denominator = G2::rand(&mut rng);
        // let witness = G1::rand(&mut rng);

        // let proof = KZGProof {
        //     numerator,
        //     denominator,
        //     witness,
        // };

        // assert!(!proof.verify(commitment));
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
