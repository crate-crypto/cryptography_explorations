use ark_bn254::Fr;
use criterion::{black_box, criterion_group, criterion_main, BatchSize, Criterion};
use cryptography_exploration::matrix_math::util::hadamard_product;
use nalgebra::DMatrix;

fn hadamard_benchmark(c: &mut Criterion) {
    let v1 = vec![
        Fr::from(1),
        Fr::from(2),
        Fr::from(3),
        Fr::from(4),
        Fr::from(5),
        Fr::from(6),
    ];
    let v2 = vec![
        Fr::from(6),
        Fr::from(5),
        Fr::from(4),
        Fr::from(3),
        Fr::from(2),
        Fr::from(1),
    ];
    let matrix = DMatrix::from_vec(1, v1.len(), v1.clone());
    let vector = DMatrix::from_vec(1, v2.len(), v2.clone());

    c.bench_function("nalgebra_benchmark", |b| {
        b.iter(|| {
            let _result = black_box(matrix.component_mul(&vector));
        })
    });

    c.bench_function("custom_benchmark", |b| {
        b.iter_batched(
            || (v1.clone(), v2.clone()),
            |(v1, v2)| black_box(hadamard_product(&v1, &v2).unwrap()),
            BatchSize::SmallInput,
        )
    });
}

criterion_group!(benches, hadamard_benchmark);
criterion_main!(benches);
