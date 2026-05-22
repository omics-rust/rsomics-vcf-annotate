use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

fn ours() -> Command {
    Command::new(env!("CARGO_BIN_EXE_rsomics-vcf-annotate"))
}

fn golden(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}

fn have(tool: &str) -> bool {
    Command::new(tool)
        .arg("--version")
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .is_ok_and(|s| s.success())
}

/// Per-variant (CHROM:POS -> ANN value or "") from a VCF.
fn ann_map(vcf: &[u8]) -> Vec<(String, String)> {
    String::from_utf8_lossy(vcf)
        .lines()
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
        .map(|l| {
            let c: Vec<&str> = l.split('\t').collect();
            let key = format!("{}:{}", c[0], c[1]);
            let ann = c
                .get(7)
                .and_then(|info| info.split(';').find(|kv| kv.starts_with("ANN=")))
                .unwrap_or("")
                .to_string();
            (key, ann)
        })
        .collect()
}

// ours' BED→INFO annotation must match `bcftools annotate -a BED -c CHROM,FROM,TO,ANN`.
#[test]
fn matches_bcftools_annotate() {
    if !have("bcftools") || !have("bgzip") || !have("tabix") {
        eprintln!("skipping: bcftools/bgzip/tabix not found");
        return;
    }
    let dir = std::env::temp_dir().join("rsomics-vcf-annotate-compat");
    let _ = std::fs::create_dir_all(&dir);

    // sort + bgzip + tabix the BED for bcftools
    let sorted = dir.join("regions.bed");
    let bed_out = std::fs::File::create(&sorted).unwrap();
    assert!(
        Command::new("sort")
            .args(["-k1,1", "-k2,2n"])
            .arg(golden("regions.bed"))
            .stdout(bed_out)
            .status()
            .unwrap()
            .success()
    );
    let gz = dir.join("regions.bed.gz");
    let gzf = std::fs::File::create(&gz).unwrap();
    assert!(
        Command::new("bgzip")
            .arg("-c")
            .arg(&sorted)
            .stdout(gzf)
            .status()
            .unwrap()
            .success()
    );
    assert!(
        Command::new("tabix")
            .args(["-fp", "bed"])
            .arg(&gz)
            .status()
            .unwrap()
            .success()
    );
    let hdr = dir.join("hdr.txt");
    std::fs::write(
        &hdr,
        "##INFO=<ID=ANN,Number=1,Type=String,Description=\"annotation\">\n",
    )
    .unwrap();

    let ours_out = ours()
        .arg(golden("small.vcf"))
        .arg("-a")
        .arg(golden("regions.bed"))
        .arg("--tag")
        .arg("ANN")
        .output()
        .unwrap();
    let bcf_out = Command::new("bcftools")
        .args(["annotate", "-a"])
        .arg(&gz)
        .args(["-c", "CHROM,FROM,TO,ANN", "-h"])
        .arg(&hdr)
        .arg(golden("small.vcf"))
        .output()
        .unwrap();
    assert!(bcf_out.status.success());

    assert_eq!(ann_map(&ours_out.stdout), ann_map(&bcf_out.stdout));
}
