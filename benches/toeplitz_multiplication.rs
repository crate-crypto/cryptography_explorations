use ark_bn254::Fr;
use ark_ff::{FftField, UniformRand};
use criterion::{criterion_group, criterion_main, Criterion};
use cryptography_exploration::toeplitz::add_zeros_to_right;
use cryptography_exploration::toeplitz::generate_toeplitz_matrix_and_vector;
use cryptography_exploration::toeplitz::ToeplitzMatrix;
use nalgebra::DMatrix;
use nalgebra::DVector;

fn toeplitz_matrix_multiplication_benchmark_4096_elements(c: &mut Criterion) {
    // Create a Toeplitz matrix and a vector for multiplication
    let mut rng: rand::prelude::ThreadRng = rand::thread_rng();
    let (toeplitz, mut vector) = generate_toeplitz_matrix_and_vector(&mut rng, 4096);
    let len = toeplitz.first_row.len();
    let mat = DMatrix::from_vec(len, len, toeplitz.get_matrix());
    let vec = DVector::from_vec(vector.to_vec());

    c.bench_function("Toeplitz Matrix Multiplication no FFT", |b| {
        b.iter(|| &mat * &vec)
    });

    let len = vector.len();
    add_zeros_to_right(&mut vector, len);

    c.bench_function("Toeplitz Matrix Multiplication FFT", |b| {
        b.iter(|| toeplitz.fast_multiply_by_vec(&vector).unwrap())
    });
}

fn toeplitz_matrix_multiplication_benchmark_8_elements(c: &mut Criterion) {
    // Create a Toeplitz matrix and a vector for multiplication
    let mut rng: rand::prelude::ThreadRng = rand::thread_rng();
    let (toeplitz, mut vector) = generate_toeplitz_matrix_and_vector(&mut rng, 8);
    let len = toeplitz.first_row.len();
    let mat = DMatrix::from_vec(len, len, toeplitz.get_matrix());
    let vec = DVector::from_vec(vector.to_vec());

    c.bench_function("Toeplitz Matrix Multiplication no FFT", |b| {
        b.iter(|| &mat * &vec)
    });

    let len = vector.len();
    add_zeros_to_right(&mut vector, len);

    c.bench_function("Toeplitz Matrix Multiplication FFT", |b| {
        b.iter(|| toeplitz.fast_multiply_by_vec(&vector).unwrap())
    });
}

criterion_group!(
    benches,
    toeplitz_matrix_multiplication_benchmark_8_elements,
    toeplitz_matrix_multiplication_benchmark_4096_elements
);
criterion_main!(benches);
