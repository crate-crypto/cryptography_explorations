use crate::{
    el_interpolation::{el_lagrange_interpolation, ElPoint},
    kzg_transcript::{KZGProof, CRS},
};
use ark_bn254::{Fr, G1Projective as G1};
use ark_ec::Group;
use ark_poly::{univariate::DensePolynomial, DenseUVPolynomial, Polynomial};
use ark_std::Zero;

#[derive(Debug, Clone, Copy)]
pub struct KZGCommitment {
    pub point: G1,
}

impl KZGCommitment {
    // TODO: MPC protocol could be implemented here
    pub fn setup(max_degree: usize) -> CRS {
        CRS::new(Fr::from(462), max_degree)
    }

    // e(q(x)1, [Z(x)]2) == e([p(x)-I(x)]1, G2)
    // numerator - I(X),
    // denominator - Z(X),
    // witness - q(x),
    // Proof: q(x) = (p(x) - I(x)) / Z(x)
    // TODO - value is here only for testing, delete it
    pub fn commit_poly(commit_to: &[ElPoint], _crs: &CRS, value: Fr) -> KZGCommitment {
        // NOTE: I checked that this commitment is generated correctly
        let commit_coeffs = el_lagrange_interpolation(commit_to);
        let commit_poly = DensePolynomial::from_coefficients_vec(commit_coeffs);
        let commitment = G1::generator() * commit_poly.evaluate(&value);

        KZGCommitment { point: commitment }
    }

    pub fn verify_poly(commitment: Self, coefficients: &[Fr], value_powers: &[G1]) -> bool {
        let mut res = G1::zero();

        for (i, coefficient) in coefficients.iter().enumerate() {
            let power = value_powers[i];
            let result = power * *coefficient;
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

    #[test]
    fn full_cycle_test() {
        let commit_to = vec![
            ElPoint::new(Fr::from(1), Fr::from(2)),
            ElPoint::new(Fr::from(2), Fr::from(3)),
            ElPoint::new(Fr::from(3), Fr::from(5)),
        ];
        let witness_to = vec![ElPoint::new(Fr::from(1), Fr::from(2))];

        let max_degree = commit_to.len();
        let value = Fr::from(1423);
        let crs = CRS::new(value, max_degree);

        // Commit to poly
        let commit_coeff = el_lagrange_interpolation(&commit_to);
        let commit_poly = DensePolynomial::from_coefficients_vec(commit_coeff.to_vec());
        let commitment = KZGCommitment::commit_poly(&commit_to, &crs, value);
        // Generate proof
        let proof = KZGProof::prove(&crs, commit_poly, &witness_to);

        // Verify polynomial
        assert!(KZGCommitment::verify_poly(
            commitment,
            &commit_coeff,
            &crs.powers_g1
        ));

        // Verify proof
        let numerator_coeffs = el_lagrange_interpolation(&witness_to);
        let numerator_poly = DensePolynomial::from_coefficients_vec(numerator_coeffs.to_vec());
        let numerator_raw = numerator_poly.evaluate(&value);
        let numerator = G1::generator() * numerator_raw;

        let zero_points: Vec<Fr> = witness_to.iter().map(|point| point.x).collect();
        let denominator_coeffs = calculate_zero_poly_coefficients(&zero_points);
        let denominator_poly = DensePolynomial::from_coefficients_vec(denominator_coeffs.to_vec());
        let denominator_raw = denominator_poly.evaluate(&value);
        let denominator = G2::generator() * denominator_raw;

        let verifier_proof = KZGProof::new(numerator, denominator, proof.witness);
        assert!(commitment.verify_proof(verifier_proof));
    }
}
