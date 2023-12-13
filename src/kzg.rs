use crate::el_interpolation::{el_lagrange_interpolation, ElPoint};
use crate::kzg_transcript::KZGParams;
use crate::kzg_transcript::KZGProof;
use crate::kzg_transcript::CRS;
use ark_bn254::{Fr, G1Projective as G1};
use ark_ec::Group;
use ark_poly::univariate::DensePolynomial;
use ark_poly::DenseUVPolynomial;
use ark_poly::Polynomial;
use ark_std::Zero;

#[derive(Debug, Clone, Copy)]
pub struct KZGCommitment {
    pub point: G1,
}

impl KZGCommitment {
    pub fn setup(value: CRS, max_degree: usize) -> KZGParams {
        KZGParams::new(max_degree, value.calc_powers(max_degree))
    }

    // e(q(x)1, [Z(x)]2) == e([p(x)-I(x)]1, G2)
    // numerator - I(X),
    // denominator - Z(X),
    // witness - q(x),
    // Proof: q(x) = (p(x) - I(x)) / Z(x)
    pub fn commit_poly(commit_to: &[ElPoint], crs: CRS) -> KZGCommitment {
        // NOTE: I checked that this commitment is generated correctly
        let commit_coeffs = el_lagrange_interpolation(commit_to);
        let commit_poly = DensePolynomial::from_coefficients_vec(commit_coeffs);
        let commitment = G1::generator() * commit_poly.evaluate(&crs.value);

        KZGCommitment { point: commitment }
    }

    pub fn verify_poly(commitment: Self, coefficients: &[G1], value_powers: &[Fr]) -> bool {
        let mut res = G1::zero();

        for (i, coefficient) in coefficients.iter().enumerate() {
            let power = value_powers[i];
            let result = *coefficient * power;
            res += result
        }

        commitment.point == res
    }

    pub fn verify_proof(&self, proof: KZGProof) -> bool {
        proof.verify(self.point)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::el_interpolation::calculate_zero_poly_coefficients;
    use ark_bn254::{Fr, G2Projective as G2};
    use rand::Rng;

    #[test]
    fn full_cycle_test() {
        let crs = CRS::new(Fr::from(13));

        let points = vec![
            ElPoint::new(Fr::from(1), Fr::from(2)),
            ElPoint::new(Fr::from(2), Fr::from(3)),
            ElPoint::new(Fr::from(3), Fr::from(5)),
        ];
        let witness_to = vec![ElPoint::new(Fr::from(1), Fr::from(2))];

        let commit_coeff = el_lagrange_interpolation(&points);
        let commit_poly = DensePolynomial::from_coefficients_vec(commit_coeff.to_vec());

        let commitment = KZGCommitment::commit_poly(&points, crs);
        let proof = KZGProof::prove(&crs, commit_poly, &witness_to);

        // Verifier part
        let numerator_coeffs = el_lagrange_interpolation(&witness_to);
        let numerator_poly = DensePolynomial::from_coefficients_vec(numerator_coeffs.to_vec());
        let numerator_raw = numerator_poly.evaluate(&crs.value);
        let numerator = G1::generator() * numerator_raw;

        let zero_points: Vec<Fr> = witness_to.iter().map(|point| point.x).collect();
        let denominator_coeffs = calculate_zero_poly_coefficients(&zero_points);
        let denominator_poly = DensePolynomial::from_coefficients_vec(denominator_coeffs.to_vec());
        let denominator_raw = denominator_poly.evaluate(&crs.value);
        let denominator = G2::generator() * denominator_raw;

        let verifier_proof = KZGProof::new(numerator, denominator, proof.witness);
        assert!(commitment.verify_proof(verifier_proof));
    }

    #[test]
    fn test_setup() {
        let mut rng: rand::prelude::ThreadRng = rand::thread_rng();
        let value = CRS::new_rand(&mut rng);
        // Generate a rand max_degree
        let max_degree = rng.gen_range(1..10);
        let params = KZGCommitment::setup(value, max_degree);
        assert!(params.max_degree == max_degree);
        assert!(params.powers.len() == max_degree + 1);
    }
}
