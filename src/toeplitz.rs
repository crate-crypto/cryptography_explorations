use ark_ff::FftField;
use ark_poly::GeneralEvaluationDomain;
use ark_poly::{domain::DomainCoeff, EvaluationDomain};
use ark_std::vec::Vec;
use std::error::Error;
use std::ops::Mul;

// TODO: implement this function
// pub fn extend_vec_to_pow_two()

// TODO: Could be implemented with as a trait with supertrait Toeplitz
#[derive(Debug)]
pub struct CirculantMatrix<F> {
    // a1, a2, ..., an (this could be also called a polynomial coefficients representation,
    // because multiplication of polynomials is the same as matrix multiplication):
    // p(x): c = fi (0,..,0, f1, ..., fd)
    vec_representation: Vec<F>, // coefficients of the polynomial what should be multiplied by another polynomial
}

impl<F: Clone> CirculantMatrix<F> {
    pub fn new(vec_representation: Vec<F>) -> Self {
        CirculantMatrix { vec_representation }
    }

    pub fn get_len(&self) -> usize {
        self.vec_representation.len()
    }

    pub fn get_matrix(&self) -> Vec<F> {
        let n = self.vec_representation.len();
        let mut matrix = Vec::with_capacity(n * n);

        for i in 0..n {
            for j in 0..n {
                let index = (j as isize - i as isize).rem_euclid(n as isize) as usize;
                matrix.push(self.vec_representation[index].clone());
            }
        }

        matrix
    }
}

impl<F: FftField + std::ops::MulAssign<F>> CirculantMatrix<F> {
    // Input (in case of Toeplitz matrix):
    // s = ([s^d-1], ..., [s], [1], [0], ..., [0]) - CRS vector with value in which the polynomial should be calculated
    // c = (0, ..., 0, f1, ..., fd) - coefficients of the polynomial which represents commitment
    //
    // Algorithm:
    // 1) y = DFT(s)
    // 2) v = DFT(c)
    // 3) u = y * v
    // 4) result = iDFT(u)
    // Link: https://alinush.github.io/2020/03/19/multiplying-a-vector-by-a-toeplitz-matrix.html
    pub fn fast_multiply_by_vec<P>(&self, vector: &[P]) -> Result<Vec<P>, Box<dyn Error>>
    where
        P: DomainCoeff<F> + std::ops::Mul<F, Output = P>,
    {
        let matrix_len = self.vec_representation.len();
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

        // y = FFT(s) - FFT(matrix)
        let y = domain.fft(&self.vec_representation);
        // v = FFT(c) - FFT(vector)
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

pub struct ToeplitzMatrix<F: FftField> {
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

    // TODO: check that it works correctly
    pub fn new_with_coeffs(coeffs: Vec<F>) -> Result<Self, Box<dyn Error>> {
        let mut zeros = vec![F::zero(); coeffs.len()];
        zeros[0] = coeffs[0];
        Ok(ToeplitzMatrix {
            first_row: coeffs,
            first_col: zeros,
        })
    }

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
    // a2N = [a0, ... , a-N-1, a0, aN-1, ..., a1] - blocking factor = 2
    // https://alinush.github.io/2020/03/19/multiplying-a-vector-by-a-toeplitz-matrix.html
    pub fn extend_to_circulant(&self) -> CirculantMatrix<F> {
        let mut result = self.first_row.clone();
        result.push(*self.first_col.first().unwrap());

        let reversed = self.first_col.iter().skip(1).rev();
        result.extend(reversed);
        CirculantMatrix::new(result)
    }

    pub fn is_power_of_two(&self) -> bool {
        self.first_col.len().is_power_of_two() && self.first_col.len().is_power_of_two()
    }

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
    // TODO: add here vector extension
    pub fn fast_multiply_by_vec(&mut self, vector: &[F]) -> Result<Vec<F>, Box<dyn Error>> {
        let circulant = self.extend_to_circulant_pow_2();

        // TODO: extend the vector
        let vec_len = vector.len();
        let mut vector: Vec<F> = if vec_len.is_power_of_two() {
            vector.to_vec()
        } else {
            let new_len = vec_len.next_power_of_two();
            let mut extended_vector = Vec::with_capacity(new_len);
            extended_vector.extend(vec![F::zero(); new_len - vec_len]);
            extended_vector.extend_from_slice(vector);
            extended_vector
        };

        let len = vector.len();
        add_zeros_to_right(&mut vector, len);

        let result = circulant.fast_multiply_by_vec(&vector).unwrap();
        Ok(result)

        // let vec_len = vector.len();
        // let circulant = toeplitz.extend_to_circulant_pow_2();
        // let circulant_len = circulant.vec_representation.len();
        // add_zeros_to_right(&mut vector, circulant_len - vec_len);
        // let vec = DVector::from_vec(vector.to_vec());
    }
}

pub fn hadamard_product<G, F, C>(v1: &[G], v2: &[F]) -> Result<Vec<C>, String>
where
    G: Mul<Output = G> + Copy,
    F: Clone + Mul<G, Output = C> + Copy,
{
    if v1.len() > v2.len() {
        return Err("The length of vector1 should be less or equal to vector2".into());
    }

    Ok(v1.iter().zip(v2).map(|(ai, bi)| *bi * *ai).collect())
}

use ark_bn254::Fr;
use ark_std::UniformRand;
use rand::Rng;
pub fn generate_toeplitz_matrix_and_vector<R: Rng>(
    rng: &mut R,
    length: usize,
) -> (ToeplitzMatrix<Fr>, Vec<Fr>) {
    // Create a sample Toeplitz matrix
    let matrix_row = (0..length).map(|_| Fr::rand(rng)).collect::<Vec<_>>();
    let mut matrix_column = Vec::new();
    matrix_column.push(matrix_row[0]);
    let extend_with = (1..length).map(|_| Fr::rand(rng)).collect::<Vec<_>>();
    matrix_column.extend_from_slice(&extend_with);
    let toeplitz = ToeplitzMatrix::new(matrix_row, matrix_column).unwrap();

    // Create a sample vector
    let vector = (0..length).map(|_| Fr::rand(rng)).collect::<Vec<_>>();

    (toeplitz, vector)
}

pub fn generate_toeplitz_matrix_and_vector_for_proptests(
    matrix_row: Vec<i64>,
    matrix_col: Vec<i64>,
    vector: Vec<i64>,
) -> (ToeplitzMatrix<Fr>, Vec<Fr>) {
    let matrix_row = matrix_row.into_iter().map(Fr::from).collect();
    let matrix_col = matrix_col.into_iter().map(Fr::from).collect();
    let vector = vector.into_iter().map(Fr::from).collect();

    let toeplitz = ToeplitzMatrix::new(matrix_row, matrix_col).expect("Invalid matrix");

    (toeplitz, vector)
}

pub fn add_zeros_to_right<T: Clone + ark_std::Zero>(vec: &mut Vec<T>, n: usize) {
    vec.resize(vec.len() + n, T::zero());
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_bn254::{Fr, G1Projective as G1};
    use ark_ec::Group;
    use ark_std::Zero;
    use nalgebra::DMatrix;

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
    #[test]
    fn test_get_circulant_matrix() {
        // Circulant matrix:
        // 1 3 2
        // 2 1 3
        // 3 2 1
        let circulant_matrix = CirculantMatrix {
            vec_representation: vec![1, 2, 3],
        };

        let matrix = circulant_matrix.get_matrix();

        let expected_matrix = vec![1, 2, 3, 3, 1, 2, 2, 3, 1];

        assert_eq!(matrix, expected_matrix);
        assert!(is_toeplitz_matrix(&matrix));
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
        let toeplitz_len = toeplitz.first_row.len();
        let mat = DMatrix::from_vec(toeplitz_len, toeplitz_len, toeplitz.get_matrix());
        let vec = DVector::from_vec(vector.to_vec());
        let expected_result_toeplitz = mat * vec;

        // Calculate expected result by direct multiplication of Circulant (extended Toeplitz to the matrix with len = n^2) matrix by vector
        let vec_len = vector.len();
        let circulant = toeplitz.extend_to_circulant_pow_2();
        let circulant_len = circulant.vec_representation.len();
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
            // NOTE: this strings are useful for debugging
            // println!(
            //     "The actual: {:?} expected_result_circulant: {:?} expected_result_toeplitz: {:?}",
            //     result, expected_result_circulant, expected_result_toeplitz
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
    fn test_hadamard_product() {
        // Integer vectors (compatible)
        let v1 = vec![1, 2, 3];
        let v2 = vec![4, 5, 6];
        let expected_result = vec![4, 10, 18];
        assert_eq!(hadamard_product(&v1, &v2).unwrap(), expected_result);

        // Integer vectors (compatible)
        let v1 = vec![1, 2, 3];
        let v2 = vec![4, 5, 6, 7, 8];
        let expected_result = vec![4, 10, 18];
        assert_eq!(hadamard_product(&v1, &v2).unwrap(), expected_result);

        // Integer vectors (incompatible)
        let v1 = vec![1, 2, 3, 4, 5];
        let v2 = vec![4, 5, 6];
        assert_eq!(
            hadamard_product(&v1, &v2).unwrap_err(),
            "The length of vector1 should be less or equal to vector2"
        );

        // Floating-point vectors
        let v1 = vec![0.5, 0.6, 0.7];
        let v2 = vec![1.5, 3.6, 0.6];
        let expected_result2 = vec![0.75, 2.16, 0.42];
        assert_eq!(hadamard_product(&v1, &v2).unwrap(), expected_result2);

        // Fr and G1 vectors
        let v1 = vec![Fr::from(1), Fr::from(2), Fr::from(3)];
        let v2 = vec![
            G1::generator() * Fr::from(1),
            G1::generator() * Fr::from(2),
            G1::generator() * Fr::from(3),
        ];
        let expected_result2: Vec<ark_ec::short_weierstrass::Projective<ark_bn254::g1::Config>> = vec![
            G1::generator() * Fr::from(1),
            G1::generator() * Fr::from(4),
            G1::generator() * Fr::from(9),
        ];
        assert_eq!(hadamard_product(&v1, &v2).unwrap(), expected_result2);
    }

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

    #[test]
    fn test_multiply_by_vector_compatible() {
        // Create a sample circulant matrix
        let matrix = CirculantMatrix::new(vec![Fr::from(1), Fr::from(2), Fr::from(3), Fr::from(4)]);

        // Create a sample vector
        let vector = vec![Fr::from(5), Fr::from(6), Fr::from(7), Fr::from(8)];

        // Perform matrix-vector multiplication
        let result = matrix.fast_multiply_by_vec(&vector).unwrap();

        assert_eq!(
            result,
            vec![Fr::from(66), Fr::from(68), Fr::from(66), Fr::from(60),]
        );
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

    #[test]
    fn test_add_zeros_to_right() {
        // Create a sample vector of Fr elements
        let mut vec = vec![Fr::from(1), Fr::from(2), Fr::from(3)];

        // Add 3 zeros to the right of the vector
        let n = 3;
        add_zeros_to_right(&mut vec, n);

        // Create the expected vector with zeros added to the right
        let expected_vec = vec![
            Fr::from(1),
            Fr::from(2),
            Fr::from(3),
            Fr::zero(),
            Fr::zero(),
            Fr::zero(),
        ];

        // Assert that the actual and expected vectors are equal
        assert_eq!(vec, expected_vec);
    }

    use nalgebra::base::DVector;

    // Generate values using random values
    #[test]
    fn test_toeplitz_matrix_vector_mut_custom_len() {
        let mut rng: rand::prelude::ThreadRng = rand::thread_rng();
        let (toeplitz, mut vector) = generate_toeplitz_matrix_and_vector(&mut rng, 11);

        test_toeplitz_matrix_vector_multiplication(toeplitz, &mut vector).unwrap()
    }

    #[test]
    fn test_toeplitz_matrix_vector_mult_4() {
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
    fn test_toeplitz_matrix_vector_mult_6() {
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

    #[test]
    fn test_toeplitz_matrix_vector_mult_8() {
        // Create a sample Toeplitz matrix
        // 1 2 3 4
        // 7 1 2 3
        // 6 7 1 2
        // 5 6 7 1
        let toeplitz = ToeplitzMatrix::new(
            vec![Fr::from(1), Fr::from(2), Fr::from(3), Fr::from(4)],
            vec![Fr::from(1), Fr::from(7), Fr::from(6), Fr::from(5)],
        )
        .unwrap();

        // Create a sample vector
        let mut vector = vec![Fr::from(8), Fr::from(9), Fr::from(10), Fr::from(11)];

        test_toeplitz_matrix_vector_multiplication(toeplitz, &mut vector).unwrap()
    }

    use proptest::collection;
    use proptest::prelude::any;
    use proptest::prelude::Arbitrary;
    use proptest::prelude::BoxedStrategy;
    use proptest::proptest;
    use proptest::strategy::Strategy;
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
                        .prop_map(
                            move |(mut matrix_row, mut matrix_col, vector)| {
                                matrix_row[0] = first_element;
                                matrix_col[0] = first_element;
                                MatVecInput::new(matrix_row, matrix_col, vector)
                            },
                        )
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
}
