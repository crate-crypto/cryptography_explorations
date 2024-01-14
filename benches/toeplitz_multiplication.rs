use criterion::{criterion_group, criterion_main, Criterion};
use cryptography_exploration::matrix_math::toeplitz::ToeplitzMatrix;
use cryptography_exploration::matrix_math::util::add_zeros_to_right;
use nalgebra::DMatrix;
use nalgebra::DVector;

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

fn toeplitz_matrix_multiplication_benchmark_4096_elements(c: &mut Criterion) {
    // Create a Toeplitz matrix and a vector for multiplication
    let mut rng: rand::prelude::ThreadRng = rand::thread_rng();
    // Note: Even with extension of the matrix and vector with zeroes to the next pow two (e.g. 2049 is extended to 4096 (in two times)) it works faster
    let (mut toeplitz, mut vector) = generate_toeplitz_matrix_and_vector(&mut rng, 4096);
    let len = toeplitz.first_row.len();
    let mat = DMatrix::from_vec(len, len, toeplitz.get_matrix());
    let vec = DVector::from_vec(vector.to_vec());

    c.bench_function("Toeplitz Matrix Multiplication no FFT (len 4096)", |b| {
        b.iter(|| &mat * &vec)
    });

    let len = vector.len();
    add_zeros_to_right(&mut vector, len);

    c.bench_function("Toeplitz Matrix Multiplication FFT (len 4096)", |b| {
        b.iter(|| toeplitz.fast_multiply_by_vec(&vector).unwrap())
    });
}

fn toeplitz_matrix_multiplication_benchmark_8_elements(c: &mut Criterion) {
    // Create a Toeplitz matrix and a vector for multiplication
    let mut rng: rand::prelude::ThreadRng = rand::thread_rng();
    let (mut toeplitz, mut vector) = generate_toeplitz_matrix_and_vector(&mut rng, 8);
    let len = toeplitz.first_row.len();
    let mat = DMatrix::from_vec(len, len, toeplitz.get_matrix());
    let vec = DVector::from_vec(vector.to_vec());

    c.bench_function("Toeplitz Matrix Multiplication no FFT (len 8)", |b| {
        b.iter(|| &mat * &vec)
    });

    let len = vector.len();
    add_zeros_to_right(&mut vector, len);

    c.bench_function("Toeplitz Matrix Multiplication FFT(len 8)", |b| {
        b.iter(|| toeplitz.fast_multiply_by_vec(&vector).unwrap())
    });
}

criterion_group!(
    benches,
    toeplitz_matrix_multiplication_benchmark_8_elements,
    toeplitz_matrix_multiplication_benchmark_4096_elements
);
criterion_main!(benches);
