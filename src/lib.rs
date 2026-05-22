use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use rsomics_common::{Result, RsomicsError};
use rsomics_intervals::{Interval, IntervalIndex, IntervalSet};

pub fn annotate_vcf(
    vcf_path: &Path,
    bed_path: &Path,
    output: &mut dyn Write,
    tag: &str,
) -> Result<u64> {
    let (index, labels) = load_annotations(bed_path)?;

    let file = File::open(vcf_path)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", vcf_path.display())))?;
    let reader = BufReader::new(file);
    let mut out = BufWriter::with_capacity(256 * 1024, output);
    let mut count: u64 = 0;

    for line in reader.lines() {
        let line = line.map_err(RsomicsError::Io)?;
        if line.starts_with('#') {
            writeln!(out, "{line}").map_err(RsomicsError::Io)?;
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 8 {
            writeln!(out, "{line}").map_err(RsomicsError::Io)?;
            continue;
        }
        let chrom = fields[0];
        let pos: u64 = fields[1].parse().unwrap_or(0);

        // 1-based VCF POS overlaps a 0-based BED region [start,end) iff the
        // 0-based position pos-1 is in it — matches bcftools annotate. Pick the
        // earliest-starting region for determinism.
        let mut label: Option<&str> = None;
        let mut best_start = u64::MAX;
        if pos > 0 {
            index.for_each_overlap(chrom, pos - 1, pos, |iv| {
                if iv.start < best_start
                    && let Some(l) = labels.get(&(iv.chrom.clone(), iv.start, iv.end))
                {
                    best_start = iv.start;
                    label = Some(l.as_str());
                }
            });
        }

        if let Some(label) = label {
            for (i, field) in fields.iter().enumerate() {
                if i > 0 {
                    write!(out, "\t").map_err(RsomicsError::Io)?;
                }
                if i == 7 {
                    if *field == "." {
                        write!(out, "{tag}={label}").map_err(RsomicsError::Io)?;
                    } else {
                        write!(out, "{field};{tag}={label}").map_err(RsomicsError::Io)?;
                    }
                } else {
                    write!(out, "{field}").map_err(RsomicsError::Io)?;
                }
            }
            writeln!(out).map_err(RsomicsError::Io)?;
        } else {
            writeln!(out, "{line}").map_err(RsomicsError::Io)?;
        }
        count += 1;
    }

    out.flush().map_err(RsomicsError::Io)?;
    Ok(count)
}

#[allow(clippy::type_complexity)]
fn load_annotations(path: &Path) -> Result<(IntervalIndex, HashMap<(String, u64, u64), String>)> {
    let file = File::open(path)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", path.display())))?;
    let reader = BufReader::new(file);
    let mut intervals: Vec<Interval> = Vec::new();
    let mut labels: HashMap<(String, u64, u64), String> = HashMap::new();

    for line in reader.lines() {
        let line = line.map_err(RsomicsError::Io)?;
        if line.starts_with('#') || line.is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 4 {
            continue;
        }
        let (chrom, label) = (f[0].to_string(), f[3].to_string());
        let start: u64 = f[1].parse().unwrap_or(0);
        let end: u64 = f[2].parse().unwrap_or(0);
        // skip degenerate (empty) BED regions — they annotate nothing
        if let Ok(iv) = Interval::new(&chrom, start, end) {
            labels.insert((chrom, start, end), label);
            intervals.push(iv);
        }
    }

    let set: IntervalSet = intervals.into_iter().collect();
    Ok((IntervalIndex::build(&set), labels))
}
