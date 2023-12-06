#[derive(Debug)]
pub struct Point {
    x: f64,
    y: f64,
}

impl Point {
    pub fn new(x: f64, y: f64) -> Self {
        Point { x, y }
    }
}

pub fn lagrange_interpolation(points: &[Point]) -> Vec<f64> {
    let mut result_polynomial = vec![0.0; points.len()];
    let mut temp_polynomial = Vec::with_capacity(points.len());

    for (i, point1) in points.iter().enumerate() {
        temp_polynomial.clear();
        temp_polynomial.push(1.0);
        let mut denominator = 1.0;

        for (j, point2) in points.iter().enumerate() {
            if i != j {
                temp_polynomial.push(0.0);
                for k in (1..temp_polynomial.len()).rev() {
                    temp_polynomial[k] -= point2.x * temp_polynomial[k - 1];
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
