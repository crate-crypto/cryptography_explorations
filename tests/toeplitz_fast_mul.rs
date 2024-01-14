use ark_bn254::Fr;
use ark_poly::domain::DomainCoeff;
use cryptography_exploration::matrix_math::toeplitz::ToeplitzMatrix;
use cryptography_exploration::matrix_math::util::add_zeros_to_right;
use nalgebra::DMatrix;
use nalgebra::DVector;
use std::error::Error;

use proptest::collection;
use proptest::prelude::*;

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
    if result_circulant[0..toeplitz_len] == expected_result_circulant.as_slice()[0..toeplitz_len]
        && result_circulant[0..toeplitz_len] == expected_result_toeplitz.as_slice()[0..toeplitz_len]
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

#[derive(Debug)]
pub struct MatVecInput {
    // The main rule: the first elements of matrix_row and matrix_col must be equal
    pub matrix_row: Vec<i64>,
    pub matrix_col: Vec<i64>,
    pub vector: Vec<i64>,
}

impl Arbitrary for MatVecInput {
    type Parameters = ();
    type Strategy = BoxedStrategy<Self>;

    fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
        collection::vec(any::<i64>(), 1..=100)
            .prop_flat_map(|vec| {
                let size = vec.len();
                let first_element = vec[0];
                (
                    collection::vec(any::<i64>(), size),
                    collection::vec(any::<i64>(), size),
                    collection::vec(any::<i64>(), size),
                )
                    .prop_map(move |(mut matrix_row, mut matrix_col, vector)| {
                        matrix_row[0] = first_element;
                        matrix_col[0] = first_element;
                        MatVecInput::new(matrix_row, matrix_col, vector)
                    })
            })
            .boxed()
    }
}

impl MatVecInput {
    pub fn new(matrix_row: Vec<i64>, matrix_col: Vec<i64>, vector: Vec<i64>) -> Self {
        MatVecInput {
            matrix_row,
            matrix_col,
            vector,
        }
    }

    pub fn get_fr_matrix(self) -> Result<(ToeplitzMatrix<Fr>, Vec<Fr>), Box<dyn Error>> {
        let matrix_row = self.matrix_row.into_iter().map(Fr::from).collect();
        let matrix_col = self.matrix_col.into_iter().map(Fr::from).collect();
        let vector = self.vector.into_iter().map(Fr::from).collect();

        let toeplitz = ToeplitzMatrix::new(matrix_row, matrix_col)?;
        Ok((toeplitz, vector))
    }
}

proptest! {
    #[test]
    fn toeplitz_matrix_vector_multiplication(input in any::<MatVecInput>()) {
        let (matrix, mut vector) = input.get_fr_matrix().unwrap();

        test_toeplitz_matrix_vector_multiplication(matrix, &mut vector).unwrap();
    }
}
