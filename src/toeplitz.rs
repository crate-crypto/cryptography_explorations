use crate::GeneralEvaluationDomain;
use ark_bn254::{Fr, G1Projective as G1};
use ark_ff::FftField;
use ark_poly::EvaluationDomain;
use ark_std::vec::Vec;
use std::error::Error;
use std::ops::Mul;

// use std::borrow::Borrow;
// use ark_ff::Field;
// use ark_std::Zero;

// Vector representation
pub struct CirculantMatrix<F: FftField> {
    // a1, a2, ..., an
    vec_representation: Vec<F>, // coefficients of p(x): c = fi (0,..,0, f1, ..., fd)
}

impl<
        F: FftField + std::ops::MulAssign<ark_ff::Fp<ark_ff::MontBackend<ark_bn254::FrConfig, 4>, 4>>,
    > CirculantMatrix<F>
{
    pub fn new(vec_representation: Vec<F>) -> Self {
        CirculantMatrix { vec_representation }
    }

    // Input:
    // s = ([s^d-1], ..., [s], [1], [0], ..., [0])

    // Algorithm:
    // 1) y = DFT (s)
    // 2) v = DFT(c)
    // 3) u = y * v * (1, V, ..., V^2d-1) //NOTE: Here I'm not sure why do we need to multiply by (1, V, ..., V^2d-1)
    // 4) h = DFT(u)
    // Link: https://alinush.github.io/2020/03/19/multiplying-a-vector-by-a-toeplitz-matrix.html
    pub fn fast_multiply_by_vector(&self, vector: &[F]) -> Result<Vec<F>, Box<dyn Error>> {
        let matrix_len = self.vec_representation.len();
        // Check that the vector is the same length as the matrix
        if matrix_len != vector.len() {
            return Err("Vector lengths are not equal".into());
        }

        let domain: GeneralEvaluationDomain<F> = GeneralEvaluationDomain::new(matrix_len).unwrap();

        // y = DFT(s) - DFT(matrix)
        let y = domain.fft(&self.vec_representation);
        // v = DFT(c) - DFT(vector)
        let v = domain.fft(vector);

        // NOTE: The next two lines are equal in terms of result, but ark function panics, because: assert_eq!(self_evals.len(), other_evals.len());
        // It probably would be better to use native ark function to multiply polynomials, especially regarding ark crate development.
        // let u = hadamard_product(&y, &v).unwrap();
        let u_evals = domain.mul_polynomials_in_evaluation_domain(&y, &v);
        let u = domain.ifft(&u_evals);

        Ok(u)
    }
}

// // TODO: Could be implemented with as a trait with supertrait Toeplitz
// #[derive(Debug)]
// pub struct CirculantMatrix<F: FftField> {
//     // TODO: the matrix could be defined as polynomial coefficients:
//     // https://core.ac.uk/download/pdf/82164134.pdf
//     // vector representation
//     elements: Vec<F>,
// }

// impl<
//         F: FftField + std::borrow::Borrow<ark_ff::Fp<ark_ff::MontBackend<ark_bn254::FrConfig, 4>, 4>>,
//     > CirculantMatrix<F>
// {
//     // Create a Circulant matrix from the input array that contains all matrix elements
//     pub fn new(matrix: &[F]) -> Result<Self, Box<dyn Error>> {
//         if is_toeplitz_matrix(matrix) {
//             Ok(CirculantMatrix {
//                 elements: matrix.to_vec(),
//             })
//         } else {
//             Err("Matrix is not Toeplitz".into())
//         }
//     }

//     // Create a Circulant matrix from the input vector that contains unique matrix elements
//     // TODO: it seems like it's a transponent matrix, while I need the original
//     // https://alinush.github.io/2020/03/19/multiplying-a-vector-by-a-toeplitz-matrix.html#multiplying-a-toeplitz-matrix-by-a-vector
//     pub fn new_from_vec(vector: &[F]) -> Self {
//         let n = vector.len();
//         let mut elements = Vec::with_capacity(n * n);

//         for i in 0..n {
//             for j in 0..n {
//                 elements.push(vector[(n + i - j) % n].clone());
//             }
//         }

//         CirculantMatrix { elements }
//     }

//     // FIXME - rewrite this code
//     // https://github.com/alisaaalehi/convolution_as_multiplication/blob/main/ConvAsMulExplained.pdf
//     pub fn multiply_by_vector(&self, v: &[F]) -> Vec<F> {
//         let m = self.elements.len();
//         let n = v.len();
//         let mut result = vec![F::zero(); m];

//         for i in 0..m {
//             for j in 0..n {
//                 result[i] += self.elements[(m + i - j) % m] * v[j];
//             }
//         }

//         result
//     }

//     pub fn multiply_toeplitz_by_vector(&self, vector: &[F]) -> Result<Vec<F>, &'static str> {
//         unimplemented!()
//     }
// }

// // s = ([s d-1], [s d-2], ..., [s], [1], [0], ..., [0])  // d of 0 at the end
// // c = (0, ..., 0, f1, ..., fn) // d of z0 at the start
// pub fn toeplitz_product_by_vec(s: &[Fr], c: &[G1]) -> Result<Vec<G1>, Box<dyn Error>> {
//     // Get the number of rows and columns in F
//     let num_elems = s.len();

//     // TODO: check that s and c are the same length

//     let domain: GeneralEvaluationDomain<Fr> = GeneralEvaluationDomain::new(num_elems).unwrap();

//     // y = DFT(s)
//     let y = domain.fft(s);
//     // v = DFT(c)
//     let v = domain.fft(c);

//     // u = y * v * (1, V, ..., V^2d-1)
//     let u = hadamard_product(&y, &v).unwrap();

//     let result = domain.ifft(&u);

//     // Ok(result)
//     unimplemented!()
// }

// c = (0, ..., 0, f1, ..., fn) // d of z0 at the start
pub fn left_pad_c() {}
// s = ([s d-1], [s d-2], ..., [s], [1], [0], ..., [0])  // d of 0 at the end
pub fn right_pad_s() {}

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

// TODO: Replace with nalgebra (?)
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
    use ark_ff::BigInt;
    use ark_std::Zero;

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
    // fn test_is_toeplitz_matrix() {
    //     // Toeplitz matrix:
    //     // 1 2 3
    //     // 4 1 2
    //     // 5 4 1
    //     let toeplitz_matrix = vec![
    //         Fr::from(1),
    //         Fr::from(2),
    //         Fr::from(3),
    //         Fr::from(4),
    //         Fr::from(1),
    //         Fr::from(2),
    //         Fr::from(5),
    //         Fr::from(4),
    //         Fr::from(1),
    //     ];
    //     assert!(is_toeplitz_matrix(&toeplitz_matrix));

    //     // Non-Toeplitz matrix:
    //     // 1 2 3
    //     // 4 5 6
    //     // 7 8 9
    //     let non_toeplitz_matrix = vec![
    //         Fr::from(1),
    //         Fr::from(2),
    //         Fr::from(3),
    //         Fr::from(4),
    //         Fr::from(5),
    //         Fr::from(6),
    //         Fr::from(7),
    //         Fr::from(8),
    //         Fr::from(9),
    //     ];
    //     assert!(!is_toeplitz_matrix(&non_toeplitz_matrix));

    //     // Circulant matrix
    //     let vector = vec![Fr::from(1), Fr::from(2), Fr::from(3)];
    //     let matrix = CirculantMatrix::new_from_vec(&vector);
    //     assert!(is_toeplitz_matrix(&matrix.elements));
    // }

    // #[test]
    // fn test_new_circulant_matrix() {
    //     // Fr vector
    //     let vector = vec![Fr::from(1), Fr::from(2), Fr::from(3)];
    //     let expected_matrix = vec![
    //         Fr::from(1),
    //         Fr::from(3),
    //         Fr::from(2),
    //         Fr::from(2),
    //         Fr::from(1),
    //         Fr::from(3),
    //         Fr::from(3),
    //         Fr::from(2),
    //         Fr::from(1),
    //     ];
    //     let matrix = CirculantMatrix::new_from_vec(&vector);
    //     assert_eq!(matrix.elements, expected_matrix);
    // }

    #[test]
    fn test_multiply_by_vector() {
        // Create a sample circulant matrix
        let matrix = CirculantMatrix::new(vec![Fr::from(1), Fr::from(2), Fr::from(3), Fr::from(4)]);

        // Create a sample vector
        let vector = vec![Fr::from(5), Fr::from(6), Fr::from(7), Fr::from(8)];

        // Perform matrix-vector multiplication
        let result = matrix.fast_multiply_by_vector(&vector).unwrap();

        assert_eq!(
            result,
            vec![Fr::from(66), Fr::from(68), Fr::from(66), Fr::from(60),]
        );
    }

    #[test]
    fn test_multiply_by_vector_unequal_sizes() {
        let matrix = CirculantMatrix::new(vec![Fr::from(1), Fr::from(2), Fr::from(3)]);
        let vector = vec![Fr::from(5), Fr::from(6), Fr::from(7), Fr::from(8)];

        assert!(matrix.fast_multiply_by_vector(&vector).is_err());
    }

    // #[test]
    // fn test_multiply_toeplitz_by_vector() {
    //     // Create a sample Toeplitz matrix
    //     let matrix = vec![Fr::from(1), Fr::from(2), Fr::from(3), Fr::from(1)];

    //     // Create a sample vector
    //     let vector = vec![Fr::from(5), Fr::from(6)];

    //     // Create a CirculantMatrix instance
    //     let circulant_matrix = CirculantMatrix {
    //         elements: matrix.clone(),
    //     };

    //     // Perform Toeplitz matrix-vector multiplication
    //     let result = circulant_matrix
    //         .multiply_toeplitz_by_vector(&matrix, &vector)
    //         .unwrap();

    //     // Expected result:
    //     assert_eq!(result, vec![Fr::from(17), Fr::from(39)]);
    // }

    // TODO: AI
    // #[test]
    // fn test_fast_multiply_by_vector() {
    //     // Create a sample ToeplitzMatrix and vector for testing
    //     let matrix = CirculantMatrix {
    //         vec_representation: vec![1, 2, 3, 4],
    //     };
    //     let vector = vec![5, 6, 7, 8];

    //     // Call the fast_multiply_by_vector function
    //     let result = matrix.fast_multiply_by_vector(&vector);

    //     // Define the expected result
    //     let expected_result = vec![70, 50, 40, 38];

    //     // Check if the result matches the expected result
    //     assert_eq!(result, Ok(expected_result));
    // }

    // #[test]
    // fn test_fast_multiply_by_vector_different_lengths() {
    //     // Create a sample ToeplitzMatrix and vector for testing
    //     let matrix = CirculantMatrix {
    //         vec_representation: vec![1, 2, 3, 4],
    //     };
    //     let vector = vec![5, 6, 7];

    //     // Call the fast_multiply_by_vector function
    //     let result = matrix.fast_multiply_by_vector(&vector);

    //     // Check if the result is an Err variant with the expected error message
    //     assert_eq!(result, Err("Vector lengths are not equal"));
    // }
}
