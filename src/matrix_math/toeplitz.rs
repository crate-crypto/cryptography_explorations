use crate::matrix_math::{add_zeros_to_right, CirculantMatrix};
use ark_ff::FftField;
use ark_poly::domain::DomainCoeff;
use ark_std::vec::Vec;
use std::error::Error;

/// ## Toeplitz Matrix
///
/// The `ToeplitzMatrix` struct represents a Toeplitz matrix, which can be constructed from vectors
/// representing its first row and first column. The first elements of both vectors must be equal.
/// The Toeplitz matrix can be extended to a circulant matrix or a power-of-two-sized Toeplitz matrix.
///
/// TODO: We define two polynomials A (t) and U(t), and we reduce the problem of
/// the matrix multiplication to the problem of multiplying the two polynomials
/// A (t) and U(t). The degree of the two polynomials is of order 0 (IV), and we
/// can use the FFT for the polynomial multiplication with cost 0 (N log N) (Aho
/// et al., 1974)
pub struct ToeplitzMatrix<F> {
    // NOTE: Such representation might be not efficient, because we need to inverse the first column
    // Matrix example:
    // a0   a-1  a-2  a-3
    // a1   a0   a-1  a-2
    // a2   a1   a0   a-1
    // a3   a2   a1   a0
    pub first_row: Vec<F>, // a0 a-1 a-2 a-3
    pub first_col: Vec<F>, // a0 a1 a2 a3
}

impl<F: FftField> ToeplitzMatrix<F> {
    pub fn new(first_row: Vec<F>, first_col: Vec<F>) -> Result<Self, Box<dyn Error>> {
        if first_row.len() != first_col.len() {
            return Err("First column and first row must have the same length".into());
        }
        if first_row.is_empty() || first_col.is_empty() {
            return Err("First column and first row must contain values".into());
        }

        let first_elem_row = first_row.first().unwrap();
        let first_elem_col = first_col.first().unwrap();
        if first_elem_row != first_elem_col {
            return Err(format!(
                "The matrix is not Toeplitz: first row element = {}, first column element = {}",
                first_elem_row, first_elem_col
            )
            .into());
        }
        Ok(ToeplitzMatrix {
            first_row,
            first_col,
        })
    }

    /// Returns the matrix dimension - the length of it's vector representation (first column).
    pub fn get_len(&self) -> usize {
        self.first_col.len()
    }

    /// Function that creates a Toeplitz matrix from a vector of coefficients, so that
    pub fn new_with_coeffs(coeffs: Vec<F>) -> Result<Self, Box<dyn Error>> {
        // TODO: add here the check that coeffs isn't empty
        let mut zeros = vec![F::zero(); coeffs.len()];
        zeros[0] = coeffs[0];
        let toeplitz = ToeplitzMatrix::new(coeffs, zeros).unwrap();
        Ok(toeplitz)
    }

    /// Function that checks that the Toeplitz matrix has the dimension of power of two
    pub fn is_power_of_two(&self) -> bool {
        self.first_col.len().is_power_of_two() && self.first_col.len().is_power_of_two()
    }

    /// Function that extends both column and row to be a valid version of the same matrix,
    /// but with dimension of power of two. This operation should be valid in terms of fast multiplication
    pub fn extend_to_toeplits_pow_2(&mut self) {
        if !self.is_power_of_two() {
            let row_pow_2_len = self.first_row.len().next_power_of_two();
            let col_pow_2_len = self.first_col.len().next_power_of_two();
            self.first_row
                .extend(vec![F::zero(); row_pow_2_len - self.first_row.len()]);
            self.first_col
                .extend(vec![F::zero(); col_pow_2_len - self.first_col.len()]);
        }
    }

    /// Returns all the matrix elements in form of vector.
    pub fn get_matrix(&self) -> Vec<F> {
        let rows = self.first_row.len();
        let cols = self.first_col.len();

        let mut full_matrix = Vec::with_capacity(rows * cols);

        for i in 0..rows {
            for j in 0..cols {
                let index = if j >= i { j - i } else { i - j };
                if j >= i {
                    full_matrix.push(self.first_row[index]);
                } else {
                    full_matrix.push(self.first_col[index]);
                }
            }
        }

        full_matrix
    }

    /// Extends the toeplitz matrix to the circulant matrix by combining the two first row and column,
    /// but first extending the given Toeplitz matrix's dimension to the power of two
    /// As a result, the dimenstion of gotten Circulant matrix would be 2*n, where n is 2^k
    pub fn extend_to_circulant_pow_2(&mut self) -> CirculantMatrix<F> {
        self.extend_to_toeplits_pow_2();

        let mut result = self.first_row.clone();
        result.push(*self.first_col.first().unwrap());

        let reversed = self.first_col.iter().skip(1).rev();
        result.extend(reversed);
        CirculantMatrix::new(result)
    }
}

impl<F: FftField> ToeplitzMatrix<F> {
    pub fn fast_multiply_by_vec<P>(&self, vector: &[P]) -> Result<Vec<P>, Box<dyn Error>>
    where
        P: DomainCoeff<F> + std::ops::Mul<F, Output = P>,
    {
        let circulant: CirculantMatrix<F> = self.into();

        // Extend vector
        let vec_len = vector.len();
        let mut vector: Vec<P> = if vec_len.is_power_of_two() {
            vector.to_vec()
        } else {
            let new_len = vec_len.next_power_of_two();
            let mut extended_vector = Vec::with_capacity(new_len);
            // NOTE: Here we extend from the left side because u1,1 in terms of polynomial multiplication represents the highest power
            extended_vector.extend_from_slice(vector);
            extended_vector.extend(vec![P::zero(); new_len - vec_len]);
            extended_vector
        };

        let len = vector.len();
        // NOTE: Only after normalizing the vector to pow of two we add zeros to the right
        add_zeros_to_right(&mut vector, len);

        let result = circulant.fast_multiply_by_vec(&vector).unwrap();
        Ok(result)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::matrix_math::is_toeplitz_matrix;
    use ark_bn254::Fr;
    use nalgebra::DMatrix;
    use nalgebra::DVector;

    #[test]
    fn test_new_with_coeffs() {
        let coeffs = vec![Fr::from(1), Fr::from(2), Fr::from(3)];
        let expected_first_row = vec![Fr::from(1), Fr::from(2), Fr::from(3)];
        let expected_first_col = vec![Fr::from(1), Fr::from(0), Fr::from(0)];

        let mut toeplitz = ToeplitzMatrix::new_with_coeffs(coeffs).unwrap();

        assert_eq!(toeplitz.first_row, expected_first_row);
        assert_eq!(toeplitz.first_col, expected_first_col);

        let circulant: CirculantMatrix<Fr> = (&toeplitz).into();
        let expected_vec_repr = vec![
            Fr::from(1),
            Fr::from(2),
            Fr::from(3),
            Fr::from(1),
            Fr::from(0),
            Fr::from(0),
        ];
        assert_eq!(circulant.first_row, expected_vec_repr);

        let circulant_pow_2 = toeplitz.extend_to_circulant_pow_2();
        let expected_vec_repr_pow_2 = vec![
            Fr::from(1),
            Fr::from(2),
            Fr::from(3),
            Fr::from(0),
            Fr::from(1),
            Fr::from(0),
            Fr::from(0),
            Fr::from(0),
        ];
        assert_eq!(circulant_pow_2.first_row, expected_vec_repr_pow_2);
    }

    #[test]
    fn test_get_toeplitz_matrix() {
        // Toeplitz matrix:
        // 1 2 3
        // 4 1 2
        // 5 4 1
        let toeplitz_matrix = ToeplitzMatrix {
            first_row: vec![Fr::from(1), Fr::from(2), Fr::from(3)],
            first_col: vec![Fr::from(1), Fr::from(4), Fr::from(5)],
        };

        let matrix = toeplitz_matrix.get_matrix();

        let expected_matrix = vec![
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

        assert_eq!(matrix, expected_matrix);
        assert!(is_toeplitz_matrix(&matrix));
    }

    // TODO: vector should be moved in another way
    fn test_toeplitz_matrix_vector_multiplication<
        T: std::clone::Clone
            + ark_std::Zero
            + std::ops::Mul<T, Output = T>
            + std::cmp::PartialEq
            + std::fmt::Debug
            + DomainCoeff<T>
            + ark_ff::FftField,
    >(
        mut toeplitz: ToeplitzMatrix<T>,
        mut vector: &mut Vec<T>,
    ) -> Result<(), Box<dyn Error>> {
        // Calculate expected result by direct multiplication of Toeplitz matrix by vector
        let toeplitz_len = toeplitz.get_len();
        let mat = DMatrix::from_vec(toeplitz_len, toeplitz_len, toeplitz.get_matrix());
        let vec = DVector::from_vec(vector.to_vec());
        let expected_result_toeplitz = mat * vec;

        // Calculate expected result by direct multiplication of Circulant (extended Toeplitz to the matrix with len = n^2) matrix by vector
        let vec_len = vector.len();
        let circulant = toeplitz.extend_to_circulant_pow_2();
        let circulant_len = circulant.first_row.len();
        add_zeros_to_right(&mut vector, circulant_len - vec_len);
        let vec = DVector::from_vec(vector.to_vec());

        let mat = DMatrix::from_vec(circulant_len, circulant_len, circulant.get_matrix());

        let expected_result_circulant = mat * &vec;

        // Calculate the actual result
        let result_circulant = circulant.fast_multiply_by_vec(&vector).unwrap();
        let result_toeplitz = toeplitz.fast_multiply_by_vec(&vector).unwrap();

        // Assert that the actual and expected result vectors are equal
        if result_circulant[0..toeplitz_len]
            == expected_result_circulant.as_slice()[0..toeplitz_len]
            && result_circulant[0..toeplitz_len]
                == expected_result_toeplitz.as_slice()[0..toeplitz_len]
            && result_circulant[0..toeplitz_len] == result_toeplitz[0..toeplitz_len]
        {
            // NOTE: this code is useful for debugging:
            // println!(
            //     "The actual (calculated with Circulant matrix) result: {:?}, The actual (calculated with Toeplitz matrix) result: {:?}, expected_result_circulant: {:?}, expected_result_toeplitz: {:?}",
            //     result_circulant, result_toeplitz, expected_result_circulant, expected_result_toeplitz
            // );
            Ok(())
        } else {
            Err(format!(
    "The actual (calculated with Circulant matrix) result: {:?}, The actual (calculated with Toeplitz matrix) result: {:?}, expected_result_circulant: {:?}, expected_result_toeplitz: {:?} vectors are not equal",
    result_circulant, result_toeplitz, expected_result_circulant, expected_result_toeplitz
).into())
        }
    }

    #[test]
    fn test_toeplitz_matrix_vector_mul_4() {
        // Create a sample Toeplitz matrix
        // 1 2
        // 3 1
        let toeplitz = ToeplitzMatrix::new(
            vec![Fr::from(1), Fr::from(2)],
            vec![Fr::from(1), Fr::from(3)],
        )
        .unwrap();

        // Create a sample vector
        let mut vector = vec![Fr::from(4), Fr::from(5)];

        test_toeplitz_matrix_vector_multiplication(toeplitz, &mut vector).unwrap()
    }

    #[test]
    fn test_toeplitz_matrix_vector_mul_6() {
        // Create a sample Toeplitz matrix
        // 1 2 3
        // 4 1 2
        // 5 4 1
        let toeplitz = ToeplitzMatrix::new(
            vec![Fr::from(1), Fr::from(2), Fr::from(3)],
            vec![Fr::from(1), Fr::from(4), Fr::from(5)],
        )
        .unwrap();

        let mut vector = vec![Fr::from(6), Fr::from(7), Fr::from(8)];

        test_toeplitz_matrix_vector_multiplication(toeplitz, &mut vector).unwrap()
    }
}
