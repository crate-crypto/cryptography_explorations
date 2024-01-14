use crate::matrix_math::hadamard_product;
use crate::matrix_math::toeplitz::ToeplitzMatrix;
use ark_ff::FftField;
use ark_poly::GeneralEvaluationDomain;
use ark_poly::{domain::DomainCoeff, EvaluationDomain};
use ark_std::vec::Vec;
use std::error::Error;

/// # Circulant Matrix
///
/// The `CirculantMatrix` struct represents a circulant matrix - a square matrix in which all rows
/// are composed of the same elements and each row is rotated one element to the right relative
/// to the preceding row. It is a particular kind of Toeplitz matrix.
///
/// The structure can be constructed from a vector of field elements (`F`), which represent the
/// matrix first column.
/// In case of Toeplitz matrix fast multiplication (with extension to circulant matrix),
/// the vec representation also could be percieved as polynomial coefficients.
pub struct CirculantMatrix<F> {
    // a1, a2, ..., an (this could be also called a polynomial coefficients representation,
    // because multiplication of polynomials is the same as matrix multiplication):
    // p(x): c = fi (0,..,0, f1, ..., fd)
    pub first_row: Vec<F>, // coefficients of the polynomial what should be multiplied by another polynomial
}

impl<F> CirculantMatrix<F> {
    /// Creates a new `CirculantMatrix` instance.
    pub fn new(first_row: Vec<F>) -> Self {
        CirculantMatrix { first_row }
    }
    /// Returns the matrix dimension - the length of it's vector representation (first column).
    pub fn get_len(&self) -> usize {
        self.first_row.len()
    }
}

impl<F: Clone> CirculantMatrix<F> {
    /// Returns all the matrix elements in form of vector.
    // TODO: fix .clone()
    pub fn get_matrix(&self) -> Vec<F> {
        let n = self.first_row.len();
        let mut matrix = Vec::with_capacity(n * n);

        for i in 0..n {
            for j in 0..n {
                let index = (j as isize - i as isize).rem_euclid(n as isize) as usize;
                matrix.push(self.first_row[index].clone());
            }
        }

        matrix
    }
}

impl<F: FftField> From<&ToeplitzMatrix<F>> for CirculantMatrix<F> {
    fn from(val: &ToeplitzMatrix<F>) -> Self {
        let mut result = val.first_row.clone();
        result.push(*val.first_col.first().unwrap());

        let reversed = val.first_col.iter().skip(1).rev();
        result.extend(reversed);
        CirculantMatrix::new(result)
    }
}

impl<F: FftField + std::ops::MulAssign<F>> CirculantMatrix<F> {
    /// Function for fast multiplication of CirculantMatrix by vector
    ///
    /// Input (in case of Toeplitz matrix):
    /// s = ([s^d-1], ..., [s], [1], [0], ..., [0]) - CRS vector with value in which the polynomial should be calculated
    /// c = (0, ..., 0, f1, ..., fd) - coefficients of the polynomial which represents commitment
    ///
    /// Algorithm:
    /// 1) y = DFT(s)
    /// 2) v = DFT(c)
    /// 3) u = y * v
    /// 4) result = iDFT(u)
    /// Link to detailed information: https://alinush.github.io/2020/03/19/multiplying-a-vector-by-a-toeplitz-matrix.html
    pub fn fast_multiply_by_vec<P>(&self, vector: &[P]) -> Result<Vec<P>, Box<dyn Error>>
    where
        P: DomainCoeff<F> + std::ops::Mul<F, Output = P>,
    {
        let matrix_len = self.first_row.len();
        // Check that the vector is the same length as the matrix
        if matrix_len > vector.len() {
            return Err(format!(
                "The matrix length is more than vector's: matrix length = {}, vector length = {}",
                matrix_len,
                vector.len()
            )
            .into());
        }

        // NOTE: If the domain is for example 12 then there would be 16 proofs (because it's the next power of 2):
        let domain: GeneralEvaluationDomain<F> = GeneralEvaluationDomain::new(matrix_len).unwrap();

        // y = FFT(s) - FFT(matrix) - evals of polynomial that are represented by matrix
        let y = domain.fft(&self.first_row);
        // v = FFT(c) - FFT(vector) - evals of polynomial that are represented by vector
        let v = domain.fft(vector);

        // NOTE: The next two lines are equal in terms of result, but ark function panics, because: assert_eq!(self_evals.len(), other_evals.len());
        // It probably would be better to use native ark function to multiply polynomials, especially regarding ark crate development.
        // Moreover, ark function does not support multiplication Fr * G1. It only supports multiplication of similar inputs.
        // TODO: it would be super to make pull request to the ark_poly with the polynomial
        // TODO: multiplication, so that it would be possible to multiply not only equal values
        // TODO: (Fr on Fr), but also Fr on G1
        // NOTE: I found out that nalgebra crate works faster that my implementation
        // let u_evals = domain.mul_polynomials_in_evaluation_domain(&y, &v);
        let u_evals = hadamard_product(&y, &v).unwrap();

        let u = domain.ifft(&u_evals);

        Ok(u)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::matrix_math::is_toeplitz_matrix;
    use ark_bn254::{Fr, G1Projective as G1};
    use ark_ec::Group;

    #[test]
    fn test_get_circulant_matrix() {
        // Circulant matrix:
        // 1 3 2
        // 2 1 3
        // 3 2 1
        let circulant_matrix = CirculantMatrix {
            first_row: vec![1, 2, 3],
        };

        let matrix = circulant_matrix.get_matrix();

        let expected_matrix = vec![1, 2, 3, 3, 1, 2, 2, 3, 1];

        assert_eq!(matrix, expected_matrix);
        assert!(is_toeplitz_matrix(&matrix));
    }

    #[test]
    fn test_toeplitz_to_circulant_conversion() {
        // Create a sample Toeplitz matrix
        let toeplitz_matrix = ToeplitzMatrix::new(
            vec![Fr::from(1), Fr::from(2), Fr::from(3)],
            vec![Fr::from(1), Fr::from(4), Fr::from(5)],
        )
        .unwrap();

        // Convert the Toeplitz matrix to a Circulant matrix
        let circulant_matrix: CirculantMatrix<Fr> = (&toeplitz_matrix).into();

        // Verify the resulting Circulant matrix
        let expected_matrix = vec![
            Fr::from(1),
            Fr::from(2),
            Fr::from(3),
            Fr::from(1),
            Fr::from(5),
            Fr::from(4),
        ];
        assert_eq!(circulant_matrix.first_row, expected_matrix);
        assert!(is_toeplitz_matrix(&circulant_matrix.get_matrix()));
    }

    #[test]
    fn test_multiply_by_vector_incompatible() {
        let matrix = CirculantMatrix::new(vec![Fr::from(1), Fr::from(2), Fr::from(3), Fr::from(4)]);
        let vector = vec![Fr::from(6), Fr::from(7), Fr::from(8)];

        assert!(matrix.fast_multiply_by_vec(&vector).is_err());
    }

    #[test]
    fn test_multiply_by_vector_group() {
        // Create a sample circulant matrix
        let matrix = CirculantMatrix::new(vec![Fr::from(1), Fr::from(2), Fr::from(3), Fr::from(4)]);

        // Create a sample vector
        let vector = vec![
            G1::generator() * Fr::from(5),
            G1::generator() * Fr::from(6),
            G1::generator() * Fr::from(7),
            G1::generator() * Fr::from(8),
        ];

        // Perform matrix-vector multiplication
        let result = matrix.fast_multiply_by_vec(&vector).unwrap();

        assert_eq!(
            result,
            vec![
                G1::generator() * Fr::from(66),
                G1::generator() * Fr::from(68),
                G1::generator() * Fr::from(66),
                G1::generator() * Fr::from(60),
            ]
        );
    }
}
