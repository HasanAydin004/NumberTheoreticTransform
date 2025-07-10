use criterion::{criterion_group, criterion_main, Criterion};
use number_theoretic_transform::{forward_ntt_recursive,ntt_based_on_omega};

fn compare_ntt_algorithms(c: &mut Criterion) {
    let input: Vec<i128> = (0..256).collect();
    let root = 17;
    let modulus = 3329;

    let mut group = c.benchmark_group("NTT Comparison");

    group.bench_function("recursive NTT radix 2", |b| {
        b.iter(|| forward_ntt_recursive(&input, root, modulus))
    });

    group.bench_function("Classic NTT", |b| {
        b.iter(|| ntt_based_on_omega(&input, root, modulus))
    });

    group.finish();
}



criterion_group!(benches, compare_ntt_algorithms);
criterion_main!(benches);
