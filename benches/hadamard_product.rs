use criterion::{black_box, criterion_group, criterion_main, BatchSize, Criterion};
use cryptography_exploration::toeplitz::hadamard_product;
use nalgebra::Matrix2;

fn nalgebra_benchmark(c: &mut Criterion) {
    let matrix = Matrix2::new(1.0, 2.0, 3.0, 4.0);
    let vector = Matrix2::new(1.0, 2.0, 3.0, 4.0);

    c.bench_function("nalgebra_component_mul", |b| {
        b.iter(|| {
            let _result = black_box(matrix.component_mul(&vector));
        })
    });
}

fn custom_benchmark(c: &mut Criterion) {
    let v3 = vec![0.5, 0.6, 0.7];
    let v4 = vec![1.5, 3.6, 0.6];

    c.bench_function("hadamard_product", |b| {
        b.iter_batched(
            || (v3.clone(), v4.clone()),
            |(v1, v2)| black_box(hadamard_product(&v1, &v2).unwrap()),
            BatchSize::SmallInput,
        )
    });
}

criterion_group!(benches, nalgebra_benchmark, custom_benchmark);
criterion_main!(benches);
