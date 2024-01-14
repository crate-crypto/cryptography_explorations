use ark_std::ops::Mul;

// TODO: this one could be more efficient using nalgebra crate
pub fn hadamard_product<G, F, C>(v1: &[G], v2: &[F]) -> Result<Vec<C>, String>
where
    G: Mul<Output = G> + Copy + Sync,
    F: Clone + Mul<G, Output = C> + Copy + Sync,
    C: Send,
{
    if v1.len() > v2.len() {
        return Err("The length of vector1 should be less or equal to vector2".into());
    }

    Ok(v1.iter().zip(v2).map(|(ai, bi)| *bi * *ai).collect())
}

// Function that is responsible for vector extension from the right side
pub fn add_zeros_to_right<T: Clone + ark_std::Zero>(vec: &mut Vec<T>, n: usize) {
    vec.resize(vec.len() + n, T::zero());
}

pub fn is_toeplitz_matrix<F: PartialEq>(matrix: &[F]) -> bool {
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

#[cfg(test)]
mod tests {
    use super::*;
    use ark_bn254::{Fr, G1Projective as G1};
    use ark_ec::Group;
    use ark_std::Zero;
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
