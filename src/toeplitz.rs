use crate::GeneralEvaluationDomain;
use ark_ff::FftField;
use ark_poly::EvaluationDomain;
use ark_std::vec::Vec;
use std::error::Error;
use std::ops::Mul;

// use std::borrow::Borrow;
// use ark_ff::Field;
// use ark_std::Zero;

// TODO: Could be implemented with as a trait with supertrait Toeplitz
#[derive(Debug)]
pub struct CirculantMatrix<F: FftField> {
    // a1, a2, ..., an (this could be also called a polynomial coefficients representation,
    // because multiplication of polynomials is the same as matrix multiplication):
    // p(x): c = fi (0,..,0, f1, ..., fd)
    vec_representation: Vec<F>, // coefficients of
}

impl<F: FftField> CirculantMatrix<F> {
    pub fn new(vec_representation: Vec<F>) -> Self {
        CirculantMatrix { vec_representation }
    }

    // c2N = [a0, a1, . . . , aN−1, a0, a−(N−1), a−(N−2), . . . , a−1]⊤
    // a0   a−1  a−2  a−3    a0   a3   a2   a1
    // a1   a0   a−1  a−2    a−3  a0   a3   a2
    // a2   a1   a0   a−1    a−2  a−3  a0   a3
    // a3   a2   a1   a0     a−1  a−2  a−3  a0
    //
    // a0   a3   a2   a1     a0   a−1  a−2  a−3
    // a−3  a0   a3   a2     a1   a0   a−1  a−2
    // a−2  a−3  a0   a3     a2   a1   a0   a−1
    // a−1  a−2  a−3  a0     a3   a2   a1   a0
    // pub fn compute_circulant_matrix(toeplitz: &ToeplitzMatrix<F>) -> Self {
    //     // 1 2 3
    //     // 5 1 2
    //     // 4 5 1
    //     let mut c2n: Vec<F> = Vec::new();
    //     c2n.extend(toeplitz.first_row.iter());
    //     println!("c2n1: {:?}", c2n);
    //     c2n.push(*toeplitz.last_row.first().unwrap());
    //     println!("c2n2: {:?}", c2n);
    //     c2n.extend(toeplitz.last_row.iter().skip(1).rev());
    //     println!("c2n1: {:?}", c2n);

    //     // c2n.extend(toeplitz.last_row.iter().rev());
    //     // println!("c2n1: {:?}", c2n);
    //     // c2n.push(*toeplitz.first_row.first().unwrap());
    //     // println!("c2n2: {:?}", c2n);
    //     // c2n.extend(toeplitz.first_row.iter().skip(1).rev());
    //     // println!("c2n3: {:?}", c2n);

    //     CirculantMatrix {
    //         vec_representation: c2n,
    //     }
    // }
}

impl<
        F: FftField + std::ops::MulAssign<ark_ff::Fp<ark_ff::MontBackend<ark_bn254::FrConfig, 4>, 4>>,
    > CirculantMatrix<F>
{
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
        let u_evals = hadamard_product(&y, &v).unwrap();
        // let u_evals = domain.mul_polynomials_in_evaluation_domain(&y, &v);
        let u = domain.ifft(&u_evals);

        Ok(u)
    }
}

pub struct ToeplitzMatrix<F: FftField> {
    // Could be presented as the first and the last row or in the form of one row and one column.
    pub first_row: Vec<F>,
    pub first_col: Vec<F>,
}

impl<F: FftField> ToeplitzMatrix<F> {
    pub fn new(first_row: Vec<F>, first_col: Vec<F>) -> Result<Self, Box<dyn Error>> {
        // TODO: rewrite
        if first_row.len() != first_col.len() {
            return Err("First column and first row must have the same length".into());
        }
        if first_row.first().unwrap() != first_col.first().unwrap() {
            return Err("The matrix is not Toeplitz".into());
        }
        Ok(ToeplitzMatrix {
            first_row, // a0 a-1 a-2 a-3
            first_col, // a0 a1 a2 a3
        })
    }

    // c2N = [a0, a1, . . . , aN −1, a0, a−(N −1), a−(N −2), . . . , a−1]⊤
    fn extend_to_circulant(&self) -> CirculantMatrix<F> {
        let mut result = self.first_row.clone();
        result.push(*self.first_col.first().unwrap());

        let reversed = self.first_col.iter().skip(1).rev();
        result.extend(reversed);

        CirculantMatrix::new(result)
    }
}

// s = ([s d-1], [s d-2], ..., [s], [1], [0], ..., [0])  // d of 0 at the end
fn add_zeros_to_right<T: Clone + ark_std::Zero>(vec: &mut Vec<T>, n: usize) {
    let zeros = vec![T::zero(); n];
    vec.extend(zeros);
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
    use ark_bn254::{Fr, G1Projective as G1};
    use ark_ec::Group;
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
        let result = matrix.fast_multiply_by_vector(&vector).unwrap();

        assert_eq!(
            result,
            vec![Fr::from(66), Fr::from(68), Fr::from(66), Fr::from(60),]
        );
    }

    #[test]
    fn test_multiply_by_vector_incompatible() {
        let matrix = CirculantMatrix::new(vec![Fr::from(1), Fr::from(2), Fr::from(3)]);
        let vector = vec![Fr::from(5), Fr::from(6), Fr::from(7), Fr::from(8)];

        assert!(matrix.fast_multiply_by_vector(&vector).is_err());
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
    use nalgebra::DMatrix;
    use nalgebra::Matrix4;
    use nalgebra::Vector4;
    #[test]
    fn test_toeplitz_matrix_vector_mult() {
        // Create a sample Toeplitz matrix
        // 1 2 3
        // 5 1 2
        // 4 5 1
        let toeplitz = ToeplitzMatrix::new(
            vec![Fr::from(1), Fr::from(2)],
            vec![Fr::from(1), Fr::from(3)],
        )
        .unwrap();

        // Create a sample vector
        let mut vector = vec![Fr::from(4), Fr::from(5)];
        let len = vector.len();
        add_zeros_to_right(&mut vector, len);

        // Perform Toeplitz matrix multiplication by vector
        let circulant = toeplitz.extend_to_circulant();
        // let circulant = CirculantMatrix::compute_circulant_matrix(&toeplitz);

        let result = circulant.fast_multiply_by_vector(&vector).unwrap();

        // a0   a−1  a−2  a−3    a0   a3   a2   a1
        // a1   a0   a−1  a−2    a−3  a0   a3   a2
        // a2   a1   a0   a−1    a−2  a−3  a0   a3
        // a3   a2   a1   a0     a−1  a−2  a−3  a0
        //
        // a0   a3   a2   a1     a0   a−1  a−2  a−3
        // a−3  a0   a3   a2     a1   a0   a−1  a−2
        // a−2  a−3  a0   a3     a2   a1   a0   a−1
        // a−1  a−2  a−3  a0     a3   a2   a1   a0

        // first_row, // a0 a-1 a-2 a-3  1 2 3
        // first_col, // a0 a1 a2 a3 1 5 4

        // Fr::from(1), Fr::from(2), Fr::from(1), Fr::from(3),
        // Fr::from(3), Fr::from(1), Fr::from(2), Fr::from(1),
        // Fr::from(1), Fr::from(3), Fr::from(1), Fr::from(2),
        // Fr::from(3), Fr::from(1), Fr::from(3), Fr::from(1),
        let mat = Matrix4::from_vec(vec![
            Fr::from(1),
            Fr::from(2),
            Fr::from(1),
            Fr::from(3),
            Fr::from(3),
            Fr::from(1),
            Fr::from(2),
            Fr::from(1),
            Fr::from(1),
            Fr::from(3),
            Fr::from(1),
            Fr::from(2),
            Fr::from(3),
            Fr::from(1),
            Fr::from(3),
            Fr::from(1),
        ]);
        let vec = Vector4::from_vec(vec![Fr::from(4), Fr::from(5), Fr::zero(), Fr::zero()]);
        // Assert that the actual and expected result vectors are equal
        let expected_result = mat * vec;
        println!("expected_result: {:?}", expected_result);
        println!("result: {:?}", result);
        let matrix = CirculantMatrix::new(vec![Fr::from(1), Fr::from(2), Fr::from(1), Fr::from(3)]);

        println!("tests: {:?}", matrix.fast_multiply_by_vector(&vector));
    }
}
