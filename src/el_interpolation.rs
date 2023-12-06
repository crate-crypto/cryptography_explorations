use ark_bn254::Fq;
use ark_std::{One, Zero};

#[derive(Debug, Clone, Copy)]
pub struct ElPoint {
    pub x: Fq,
    pub y: Fq,
}

impl ElPoint {
    pub fn new(x: Fq, y: Fq) -> Self {
        ElPoint { x, y }
    }
}

pub fn el_lagrange_interpolation(points: &[ElPoint]) -> Vec<Fq> {
    let mut result_polynomial = vec![Fq::zero(); points.len()];
    let mut temp_polynomial = Vec::with_capacity(points.len());

    for (i, point1) in points.iter().enumerate() {
        temp_polynomial.clear();
        temp_polynomial.push(Fq::one());
        let mut denominator = Fq::one();

        for (j, point2) in points.iter().enumerate() {
            if i != j {
                temp_polynomial.push(Fq::zero());
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

    let reversed_coefficients: Vec<_> = result_polynomial.iter().copied().rev().collect();

    reversed_coefficients
}
