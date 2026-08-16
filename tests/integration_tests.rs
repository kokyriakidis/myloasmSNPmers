use assert_cmd::prelude::*;
use std::fs;
use std::path::Path;
use std::path::PathBuf;
use std::process::Command;
use tempfile::TempDir;

// SNPmer-only fork: these tests exercise the single remaining function of the
// binary -- run k-mer counting + SNPmer detection and write snpmers.tsv. The
// original assembly / checkpoint-resume tests were removed along with the
// downstream pipeline.

/// Run the binary on `test_input` into a fresh temp output dir and return the
/// path to the produced snpmers.tsv (asserting the run succeeded and the file
/// exists).
fn run_snpmers(
    test_input: &Path,
    extra_args: &[&str],
) -> Result<(TempDir, PathBuf), Box<dyn std::error::Error>> {
    let temp_dir = TempDir::new()?;
    let output_dir = temp_dir.path().join("output");

    let mut cmd = Command::cargo_bin("myloasm")?;
    cmd.arg(test_input.to_str().unwrap())
        .arg("-o")
        .arg(output_dir.to_str().unwrap())
        .arg("-t")
        .arg("4");
    for a in extra_args {
        cmd.arg(a);
    }
    cmd.assert().success();

    let snpmers = output_dir.join("snpmers.tsv");
    assert!(
        snpmers.exists(),
        "snpmers.tsv was not produced at {:?}",
        snpmers
    );
    Ok((temp_dir, snpmers))
}

/// Parse snpmers.tsv (checking the header, then skipping it) into rows.
fn read_snpmers(path: &Path) -> Vec<Vec<String>> {
    let content = fs::read_to_string(path).expect("read snpmers.tsv");
    let mut lines = content.lines();
    let header = lines.next().expect("snpmers.tsv is empty (no header)");
    assert_eq!(
        header,
        "split_kmer\tmid_pos\tallele0_base\tallele1_base\tallele0_kmer\tallele1_kmer\tallele0_count\tallele1_count",
        "unexpected snpmers.tsv header"
    );
    lines
        .filter(|l| !l.is_empty())
        .map(|l| l.split('\t').map(|s| s.to_string()).collect())
        .collect()
}

/// Every SNPmer must be a pair of k-mers that differ at exactly one position,
/// and that position must be the reported middle position; the middle base of
/// each allele k-mer must match the reported allele base. An empty set is valid
/// (e.g. a small homozygous input), so this only validates the rows present.
fn assert_snpmers_wellformed(rows: &[Vec<String>]) {
    for r in rows {
        assert_eq!(r.len(), 8, "unexpected column count in row: {:?}", r);
        let mid_pos: usize = r[1].parse().expect("mid_pos");
        let a0 = &r[4];
        let a1 = &r[5];
        assert_eq!(a0.len(), a1.len(), "allele k-mers differ in length");

        let bytes0 = a0.as_bytes();
        let bytes1 = a1.as_bytes();
        let mut diffs = 0usize;
        let mut diff_at_mid = false;
        for i in 0..bytes0.len() {
            if bytes0[i] != bytes1[i] {
                diffs += 1;
                if i == mid_pos {
                    diff_at_mid = true;
                }
            }
        }
        assert_eq!(diffs, 1, "SNPmer alleles must differ at exactly one base: {:?}", r);
        assert!(diff_at_mid, "the single difference must be at mid_pos: {:?}", r);

        // Reported allele bases must match the middle base of each k-mer.
        assert_eq!(
            &a0[mid_pos..mid_pos + 1],
            r[2],
            "allele0 mid base mismatch: {:?}",
            r
        );
        assert_eq!(
            &a1[mid_pos..mid_pos + 1],
            r[3],
            "allele1 mid base mismatch: {:?}",
            r
        );

        // Counts must be positive integers.
        let c0: u64 = r[6].parse().expect("allele0_count");
        let c1: u64 = r[7].parse().expect("allele1_count");
        assert!(c0 > 0 && c1 > 0, "allele counts must be positive: {:?}", r);
    }
}

#[test]
fn test_missing_input() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("myloasm")?;
    cmd.assert().failure();
    Ok(())
}

#[test]
fn test_snpmers_produced_2kb_plasmid() -> Result<(), Box<dyn std::error::Error>> {
    let (_tmp, snpmers) = run_snpmers(Path::new("tests/reads/2kb_plas.fq"), &[])?;
    let rows = read_snpmers(&snpmers);
    assert_snpmers_wellformed(&rows);
    Ok(())
}

#[test]
fn test_snpmers_produced_48kb_plasmid() -> Result<(), Box<dyn std::error::Error>> {
    let (_tmp, snpmers) = run_snpmers(Path::new("tests/reads/40kb_plas.fq"), &[])?;
    let rows = read_snpmers(&snpmers);
    assert!(!rows.is_empty(), "expected SNPmers from the 48kb plasmid");
    assert_snpmers_wellformed(&rows);
    Ok(())
}

#[test]
fn test_snpmers_produced_48kb_plasmid_fasta() -> Result<(), Box<dyn std::error::Error>> {
    let (_tmp, snpmers) = run_snpmers(Path::new("tests/reads/40kb_plas.fa"), &[])?;
    let rows = read_snpmers(&snpmers);
    assert_snpmers_wellformed(&rows);
    Ok(())
}

#[test]
fn test_snpmers_deterministic() -> Result<(), Box<dyn std::error::Error>> {
    // Two runs on the same input must produce identical snpmers.tsv.
    let (_tmp1, s1) = run_snpmers(Path::new("tests/reads/40kb_plas.fq"), &[])?;
    let (_tmp2, s2) = run_snpmers(Path::new("tests/reads/40kb_plas.fq"), &[])?;
    let c1 = fs::read_to_string(&s1)?;
    let c2 = fs::read_to_string(&s2)?;
    assert_eq!(c1, c2, "snpmers.tsv differed between two runs on the same input");
    Ok(())
}
