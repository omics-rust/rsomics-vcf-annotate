use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;
use std::path::PathBuf;
use std::process::Command;

fn bench_vcf_annotate(c: &mut Criterion) {
    let bin = env!("CARGO_BIN_EXE_rsomics-vcf-annotate");
    let manifest = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let vcf = manifest.join("tests/golden/small.vcf");
    let ann = manifest.join("tests/golden/regions.bed");
    c.bench_function("rsomics-vcf-annotate golden", |b| {
        b.iter(|| {
            let out = Command::new(black_box(bin))
                .args([
                    vcf.to_str().unwrap(),
                    "-a",
                    ann.to_str().unwrap(),
                    "-o",
                    "/dev/null",
                ])
                .output()
                .unwrap();
            assert!(out.status.success());
        });
    });
}

criterion_group!(benches, bench_vcf_annotate);
criterion_main!(benches);
