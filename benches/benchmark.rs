use criterion::BatchSize;
use criterion::{criterion_group, criterion_main, Criterion};
use nalgebra::{Matrix2, Vector2};
use std::error::Error;
use std::ops::Mul;

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

fn nalgebra_benchmark(c: &mut Criterion) {
    let matrix = Matrix2::new(1.0, 2.0, 3.0, 4.0);
    let vector = Matrix2::new(1.0, 2.0, 3.0, 4.0);

    c.bench_function("nalgebra_component_mul", |b| {
        b.iter(|| {
            let _result = matrix.component_mul(&vector);
        })
    });
}

fn custom_benchmark(c: &mut Criterion) {
    let v3 = vec![0.5, 0.6, 0.7];
    let v4 = vec![1.5, 3.6, 0.6];

    c.bench_function("hadamard_product", |b| {
        b.iter_batched(
            || (v3.clone(), v4.clone()),
            |(v1, v2)| hadamard_product(&v1, &v2).unwrap(),
            BatchSize::SmallInput,
        )
    });
}

// fn custom_benchmark(c: &mut Criterion) {
//     let matrix = Matrix2::new(1.0, 2.0, 3.0, 4.0);
//     let vector = Vector2::new(5.0, 6.0);

//     c.bench_function("custom_hadamard_product", |b| {
//         b.iter(|| {
//             let _result = hadamard_product(&matrix, &vector);
//         })
//     });
// }

criterion_group!(benches, nalgebra_benchmark, custom_benchmark);
criterion_main!(benches);
