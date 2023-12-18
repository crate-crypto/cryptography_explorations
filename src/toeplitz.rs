use crate::GeneralEvaluationDomain;
use ark_bn254::{Fr, G1Projective as G1};
use ark_ff::FftField;
use ark_poly::EvaluationDomain;
use ark_std::vec::Vec;
use std::error::Error;
use std::ops::Mul;
// use nalgebra::base::Matrix;
// use std::borrow::Borrow;
// use ark_ff::Field;
// use ark_std::Zero;

// The generic Matrix type has four type parameters:
// T: for the matrix components scalar type.
// R: for the matrix number of rows.
// C: for the matrix number of columns.
// S: for the matrix data storage, i.e., the buffer that actually contains the matrix components.

pub struct CirculantMatrix<F: FftField> {
    // Should the whole matrix be stored here (?)
    elements: Vec<F>,
}

impl<
        F: FftField + std::borrow::Borrow<ark_ff::Fp<ark_ff::MontBackend<ark_bn254::FrConfig, 4>, 4>>,
    > CirculantMatrix<F>
{
    // Create a Toeplitz matrix from the input array that contains all matrix elements
    pub fn new(matrix: &[F]) -> Result<Self, Box<dyn Error>> {
        if is_toeplitz_matrix(matrix) {
            Ok(CirculantMatrix {
                elements: matrix.to_vec(),
            })
        } else {
            Err("Matrix is not Toeplitz".into())
        }
    }

    // FIXME - rewrite this code
    pub fn multiply_by_vector(&self, v: &[G1]) -> Vec<G1> {
        let _m = self.elements.len();
        let _n = v.len();

        unimplemented!()
    }
}

// https://alinush.github.io/2020/03/19/multiplying-a-vector-by-a-toeplitz-matrix.html
// s = ([s d-1], [s d-2], ..., [s], [1], [0], ..., [0])  // d of 0 at the end
// c = (0, ..., 0, f1, ..., fn) // d of z0 at the start
pub fn toeplitz_product(s: &[Fr], c: &[G1]) -> Result<Vec<G1>, Box<dyn Error>> {
    // Get the number of rows and columns in F
    let num_elems = s.len();
    let domain: GeneralEvaluationDomain<Fr> = GeneralEvaluationDomain::new(num_elems).unwrap();

    // y = DFT(s)
    let y = domain.fft(s);
    // v = DFT(c)
    let v = domain.fft(c);

    // let u = y*v*(1, q^1, ..., q^2d-1)
    hadamard_product(&y, &v)
}

fn is_toeplitz_matrix<F: PartialEq>(matrix: &[F]) -> bool {
    let n = (matrix.len() as f64).sqrt() as usize;

    for i in 0..n - 1 {
        for j in 0..n - 1 {
            if matrix[i * n + j] != matrix[(i + 1) * n + j + 1] {
                return false;
            }
        }
    }

    true
}

pub fn poly_to_toeplitz(_coeffs: &[Fr]) {}

pub fn hadamard_product<G, F, C>(v1: &[G], v2: &[F]) -> Result<Vec<C>, Box<dyn Error>>
where
    G: Mul<Output = G> + Copy,
    F: Clone + Mul<G, Output = C>,
    Vec<G>: FromIterator<G>,
{
    if v1.len() != v2.len() {
        return Err("Vector lengths must match".into());
    };

    Ok(v1
        .iter()
        .zip(v2.iter())
        .map(|(ai, bi)| bi.clone() * *ai)
        .collect())
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_ec::Group;
    // use ark_std::UniformRand;
    // use ark_std::Zero;

    #[test]
    fn test_hadamard_product() {
        // Integer vectors
        let v1 = vec![1, 2, 3];
        let v2 = vec![4, 5, 6];
        let expected_result = vec![4, 10, 18];
        assert_eq!(hadamard_product(&v1, &v2).unwrap(), expected_result);

        // Floating-point vectors
        let v3 = vec![0.5, 0.6, 0.7];
        let v4 = vec![1.5, 3.6, 0.6];
        let expected_result2 = vec![0.75, 2.16, 0.42];
        assert_eq!(hadamard_product(&v3, &v4).unwrap(), expected_result2);

        // Fr and G1 vectors
        let v5 = vec![Fr::from(1), Fr::from(2), Fr::from(3)];
        let v6 = vec![
            G1::generator() * Fr::from(1),
            G1::generator() * Fr::from(2),
            G1::generator() * Fr::from(3),
        ];
        let expected_result2: Vec<ark_ec::short_weierstrass::Projective<ark_bn254::g1::Config>> = vec![
            G1::generator() * Fr::from(1),
            G1::generator() * Fr::from(4),
            G1::generator() * Fr::from(9),
        ];
        assert_eq!(hadamard_product(&v5, &v6).unwrap(), expected_result2);
    }

    // #[test]
    // fn test_multiply_by_vector() {
    //     let mut rng = ark_std::test_rng();
    //     let elements = vec![Fr::rand(&mut rng), Fr::rand(&mut rng)];
    //     let circulant_matrix = CirculantMatrix { elements };

    //     let vector = vec![G1::zero(), G1::zero()];

    //     let result = circulant_matrix.multiply_by_vector(&vector);

    //     let expected_result = vec![G1::zero(), G1::zero()];
    //     assert_eq!(result, expected_result);
    // }

    #[test]
    fn test_is_toeplitz_matrix() {
        // Toeplitz matrix:
        // 1 2 3
        // 4 1 2
        // 5 4 1
        let toeplitz_matrix = vec![
            Fr::from(1),
            Fr::from(2),
            Fr::from(3),
            Fr::from(4),
            Fr::from(1),
            Fr::from(2),
            Fr::from(5),
            Fr::from(4),
            Fr::from(1),
        ];
        assert!(is_toeplitz_matrix(&toeplitz_matrix));

        // Non-Toeplitz matrix:
        // 1 2 3
        // 4 5 6
        // 7 8 9
        let non_toeplitz_matrix = vec![
            Fr::from(1),
            Fr::from(2),
            Fr::from(3),
            Fr::from(4),
            Fr::from(5),
            Fr::from(6),
            Fr::from(7),
            Fr::from(8),
            Fr::from(9),
        ];
        assert!(!is_toeplitz_matrix(&non_toeplitz_matrix));
    }
}
