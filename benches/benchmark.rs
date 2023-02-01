use criterion::{criterion_group, criterion_main, Criterion};

//fn create_metrics(samples: Vec<&[u8]>, gq: i32, freq: f32, population: &str) {
    //let _ = pqc::Metrics::new(samples, gq, freq, population);
//}

//fn add_counts(mut metrics: pqc::Metrics) {
    //let a = [1,1,1,1,1,1,1];
    //metrics.add_raw(0, &a);
//}

fn check_is_variant(gts: &Vec<&[i32]>) {
    for gt in gts {
        pqc::is_variant(gt);
    }
}

fn bench_create_data(c: &mut Criterion) {
    let mut group = c.benchmark_group("create_data");
    //for i in [1000, 10000].iter() {
    //    group.bench_with_input(BenchmarkId::new("create_metrics", i), i, 
    //        |b, i| b.iter(|| create_metrics(vec![b"test"; *i], 20, 0.001, "topmed")));
    //    group.bench_with_input(BenchmarkId::new("add_counts", i), i, 
    //        |b, i| b.iter(|| add_counts(pqc::Metrics::new(vec![b"test"; *i], 20, 0.001, "topmed"))));
    //}

    let mut genotypes: Vec<&[i32]> = vec![];
    for _ in 0..5000 {
        genotypes.push(&[2,2]);
    }
    for _ in 0..5000 {
        genotypes.push(&[2,4]);
    }
    for _ in 0..100 {
        genotypes.push(&[0,2]);
    }
    for _ in 0..1000 {
        genotypes.push(&[0,4]);
    }
    group.bench_function("check_is_variant", |b| b.iter(|| check_is_variant(&genotypes)));
    group.finish();
}

criterion_group!(benches, bench_create_data);
criterion_main!(benches);



