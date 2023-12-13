use ark_bn254::Fr;
use ark_ff::Field;
use ark_poly::univariate::DensePolynomial;
use ark_poly::DenseUVPolynomial;
use ark_std::{One, Zero};
use std::ops::Mul;

// TODO: For now it's just simpler to use it in this way, but it should be fixed in the future
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ElPoint {
    pub x: Fr,
    pub y: Fr,
}

impl ElPoint {
    pub fn new(x: Fr, y: Fr) -> Self {
        ElPoint { x, y }
    }
}
pub fn el_lagrange_interpolation(points: &[ElPoint]) -> Vec<ark_bn254::Fr> {
    let mut result_polynomial = vec![ark_bn254::Fr::zero(); points.len()];
    let mut temp_polynomial = Vec::with_capacity(points.len());

    for (i, point1) in points.iter().enumerate() {
        temp_polynomial.clear();
        temp_polynomial.push(ark_bn254::Fr::one());
        let mut denominator = ark_bn254::Fr::one();

        for (j, point2) in points.iter().enumerate() {
            if i != j {
                temp_polynomial.push(ark_bn254::Fr::zero());
                let temp_polynomial_clone = temp_polynomial.clone();
                for k in (1..temp_polynomial.len()).rev() {
                    temp_polynomial[k] -= point2.x * temp_polynomial_clone[k - 1];
                }
                denominator *= point1.x - point2.x;
            }
        }

        let multiplier = point1.y / denominator;

        for (result_coefficient, temp_coefficient) in result_polynomial
            .iter_mut()
            .zip(temp_polynomial.iter().copied())
        {
            *result_coefficient += multiplier * temp_coefficient;
        }
    }

    result_polynomial.iter().copied().rev().collect()
}

pub fn calculate_zero_poly_coefficients(roots: &[Fr]) -> Vec<Fr> {
    let mut coefficients = Vec::new();

    for &root in roots {
        // For each root, create a factor (x - root)
        let factor = DensePolynomial::from_coefficients_vec(vec![Fr::one(), -root]);

        // Multiply the current coefficients by the new factor
        coefficients = if coefficients.is_empty() {
            factor.coeffs
        } else {
            DensePolynomial::from_coefficients_vec(coefficients)
                .mul(&factor)
                .coeffs
        };
    }

    // The result is the coefficients of the zero polynomial
    coefficients.iter().copied().rev().collect()
}

pub fn dense_to_sparse<F: Zero + PartialEq + Clone>(coefficients: &[F]) -> Vec<(usize, F)> {
    let mut sparse_coeffs = Vec::new();

    for (degree, &ref coefficient) in coefficients.iter().enumerate() {
        if *coefficient != F::zero() {
            sparse_coeffs.push((degree, coefficient.clone()));
        }
    }

    sparse_coeffs
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_poly::Polynomial;

    #[test]
    fn test_el_lagrange_interpolation_empty() {
        // Define the input points and expected coefficients
        let points = vec![];

        // Call the function to get the actual coefficients
        let coefficients = el_lagrange_interpolation(&points);
        assert!(coefficients.is_empty())
    }

    #[test]
    fn test_el_lagrange_interpolation() {
        // Define the input points and expected coefficients
        let points = vec![
            ElPoint::new(Fr::from(1), Fr::from(2)),
            ElPoint::new(Fr::from(2), Fr::from(3)),
            ElPoint::new(Fr::from(3), Fr::from(5)),
            ElPoint::new(Fr::from(6), Fr::from(46)),
            ElPoint::new(Fr::from(8), Fr::from(16)),
            ElPoint::new(Fr::from(10), Fr::from(64)),
        ];

        // Call the function to get the actual coefficients
        let coefficients = el_lagrange_interpolation(&points);
        // Calculated polynomial
        let poly = DensePolynomial::from_coefficients_vec(coefficients);

        for point in points {
            assert_eq!(point.y, poly.evaluate(&point.x));
        }
    }

    #[test]
    fn test_calculate_zero_poly_coefficients1() {
        // Empty roots
        let roots: Vec<Fr> = vec![];
        let expected_coefficients: Vec<Fr> = vec![];
        let coefficients = calculate_zero_poly_coefficients(&roots);
        assert_eq!(coefficients, expected_coefficients);

        // Single root
        let roots = vec![Fr::from(1)];
        let expected_coefficients = vec![Fr::from(-1), Fr::from(1)];
        let coefficients = calculate_zero_poly_coefficients(&roots);
        assert_eq!(coefficients, expected_coefficients);

        // Multiple roots
        let roots = vec![Fr::from(1), Fr::from(2), Fr::from(3)];
        let expected_coefficients = vec![Fr::from(-6), Fr::from(11), Fr::from(-6), Fr::from(1)];
        let coefficients = calculate_zero_poly_coefficients(&roots);
        assert_eq!(coefficients, expected_coefficients);
    }
}
