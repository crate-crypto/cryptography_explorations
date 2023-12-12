use ark_bn254::Fr;
use ark_poly::univariate::DensePolynomial;
use ark_poly::DenseUVPolynomial;
use ark_std::{One, Zero};
use std::ops::Mul;

// TODO: For now it's just simpler to use it in this way, but it should be fixed in the future
#[derive(Debug, Clone, Copy)]
pub struct ElPoint {
    pub x: Fr,
    pub y: Fr,
}

impl ElPoint {
    pub fn new(x: Fr, y: Fr) -> Self {
        ElPoint { x, y }
    }
}

pub fn el_lagrange_interpolation(points: &[ElPoint]) -> Vec<Fr> {
    let mut result_polynomial = vec![Fr::zero(); points.len()];
    let mut temp_polynomial = Vec::with_capacity(points.len());

    for (i, point1) in points.iter().enumerate() {
        temp_polynomial.clear();
        temp_polynomial.push(Fr::one());
        let mut denominator = Fr::one();

        for (j, point2) in points.iter().enumerate() {
            if i != j {
                temp_polynomial.push(Fr::zero());
                for k in (1..temp_polynomial.len()).rev() {
                    let clone = &temp_polynomial[k - 1].clone();
                    temp_polynomial[k] -= point2.x * clone;
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

#[cfg(test)]
mod tests {
    use super::*;
    use ark_test_curves::BigInt;

    #[test]
    fn test_el_lagrange_interpolation() {
        // Define the input points and expected coefficients
        let points = vec![
            ElPoint::new(Fr::from(1), Fr::from(2)),
            ElPoint::new(Fr::from(2), Fr::from(3)),
            ElPoint::new(Fr::from(3), Fr::from(5)),
        ];

        // TODO: Check that everything is correct here
        let expected_coefficients = vec![
            Fr::from(2),
            // -0.5
            Fr::from(BigInt([
                11669102379873075200,
                10671829228508198984,
                15863968012492123182,
                1743499133401485332,
            ])),
            // 0.5
            Fr::from(BigInt([
                11669102379873075201,
                10671829228508198984,
                15863968012492123182,
                1743499133401485332,
            ])),
        ];

        // Call the function to get the actual coefficients
        let actual_coefficients = el_lagrange_interpolation(&points);

        // Assert that the actual coefficients match the expected coefficients
        assert_eq!(actual_coefficients, expected_coefficients);
    }
    #[test]
    fn test_calculate_zero_poly_coefficients() {
        let numbers = vec![Fr::from(1), Fr::from(2), Fr::from(3)];
        let expected_coefficients = vec![Fr::from(-6), Fr::from(11), Fr::from(-6), Fr::from(1)];

        let coefficients = calculate_zero_poly_coefficients(&numbers);

        assert_eq!(coefficients, expected_coefficients);
    }
}
