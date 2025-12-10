use assert_cmd::Command;
use predicates::prelude::*;
use std::fs::{self, File};
use std::io::Write;
use std::path::Path; // wrapper providing write_stdin()

use flate2::write::GzEncoder;
use flate2::Compression;

use bio::alphabets::dna;

/// Minimal FASTA with two records; second has a description.
const FASTA: &str = "\
>readA some desc
ACGTACGTAC
>readB
GGGCCCaaattt
";

/// TSV mapping orientations: readA -> '+', readB -> '-'
const TSV: &str = "\
ReadName\torientation
readA\t+
readB\t-
";

/// TSV with an invalid orientation token
const BAD_TSV: &str = "\
ReadName\torientation
readA\twut
";

fn write(path: &Path, s: &str) {
    fs::write(path, s).unwrap();
}

fn write_gz(path: &Path, s: &str) {
    let f = File::create(path).unwrap();
    let mut gz = GzEncoder::new(f, Compression::default());
    gz.write_all(s.as_bytes()).unwrap();
    gz.finish().unwrap();
}

fn run_ok(cmd: &mut Command) -> String {
    let out = cmd.assert().success().get_output().stdout.clone();
    String::from_utf8(out).unwrap()
}

#[test]
fn keeps_plus_flips_minus_with_suffix_and_wrap() {
    let td = tempfile::tempdir().unwrap();
    let fasta_p = td.path().join("in.fa");
    let tsv_p = td.path().join("map.tsv");

    write(&fasta_p, FASTA);
    write(&tsv_p, TSV);

    // Use a long sequence to test 60-col wrapping
    let long_seq = "A".repeat(130);
    fs::write(
        &fasta_p,
        format!(">readA some desc\n{}\n>readB\n{}\n", long_seq, long_seq),
    )
    .unwrap();

    // target + ; readB is '-' so it should flip and get suffix
    let mut cmd = Command::cargo_bin("restrand-fasta").unwrap();
    let out = run_ok(cmd.args([
        "-f",
        fasta_p.to_str().unwrap(),
        "-t",
        tsv_p.to_str().unwrap(),
        "--target-orientation",
        "+",
        "--flipped-suffix",
        "/rc",
    ]));

    // Record 1 header preserved (not flipped)
    assert!(out.lines().next().unwrap().starts_with(">readA some desc"));

    // Record 1 wrapping: 130 nt -> 60 + 60 + 10
    let lines: Vec<&str> = out.lines().collect();
    // headers at 0 and 3 (because 1..=3 lines are seq)
    assert_eq!(lines[1].len(), 60);
    assert_eq!(lines[2].len(), 60);
    assert_eq!(lines[3].len(), 10);

    // Record 2 header has suffix
    let rec2_header_idx = lines.iter().position(|l| l.starts_with(">readB")).unwrap();
    assert!(lines[rec2_header_idx].contains("/rc"));

    // Record 2 sequence is reverse-complement of long_seq
    let expected_rc = dna::revcomp(long_seq.as_bytes());
    let rec2_seq: String = lines[(rec2_header_idx + 1)..]
        .iter()
        .take_while(|l| !l.starts_with('>'))
        .copied()
        .collect();
    assert_eq!(rec2_seq.as_bytes(), expected_rc.as_slice());
}

#[test]
fn drop_missing_false_passes_through() {
    let td = tempfile::tempdir().unwrap();
    let fasta_p = td.path().join("in.fa");
    let tsv_p = td.path().join("map.tsv");
    write(&fasta_p, FASTA);
    // Only map readA; readB missing
    write(&tsv_p, "ReadName\torientation\nreadA\t+\n");

    let mut cmd = Command::cargo_bin("restrand-fasta").unwrap();
    let out = run_ok(cmd.args([
        "-f",
        fasta_p.to_str().unwrap(),
        "-t",
        tsv_p.to_str().unwrap(),
    ]));

    // Both records should be present (readB passed through)
    assert!(out.contains(">readA some desc"));
    assert!(out.contains(">readB"));
}

#[test]
fn drop_missing_true_drops() {
    let td = tempfile::tempdir().unwrap();
    let fasta_p = td.path().join("in.fa");
    let tsv_p = td.path().join("map.tsv");
    write(&fasta_p, FASTA);
    // Only map readA; readB missing
    write(&tsv_p, "ReadName\torientation\nreadA\t+\n");

    let mut cmd = Command::cargo_bin("restrand-fasta").unwrap();
    let out = run_ok(cmd.args([
        "-f",
        fasta_p.to_str().unwrap(),
        "-t",
        tsv_p.to_str().unwrap(),
        "--drop-missing",
    ]));

    // readB should be gone
    assert!(out.contains(">readA some desc"));
    assert!(!out.contains(">readB"));
}

#[test]
fn gz_fasta_and_gz_table_work() {
    let td = tempfile::tempdir().unwrap();
    let fasta_gz = td.path().join("in.fa.gz");
    let tsv_gz = td.path().join("map.tsv.gz");

    write_gz(&fasta_gz, FASTA);
    write_gz(&tsv_gz, TSV);

    let mut cmd = Command::cargo_bin("restrand-fasta").unwrap();
    let out = run_ok(cmd.args([
        "-f",
        fasta_gz.to_str().unwrap(),
        "-t",
        tsv_gz.to_str().unwrap(),
    ]));
    // basic sanity: both headers present
    assert!(out.contains(">readA some desc"));
    assert!(out.contains(">readB"));
}

#[test]
fn stdin_stdout_mode() {
    let td = tempfile::tempdir().unwrap();
    let tsv_p = td.path().join("map.tsv");
    write(&tsv_p, TSV);

    // feed FASTA on stdin
    let mut cmd = Command::cargo_bin("restrand-fasta").unwrap();
    cmd.args(["-f", "-", "-t", tsv_p.to_str().unwrap()]);
    cmd.write_stdin(FASTA);

    let out = run_ok(&mut cmd);
    assert!(out.contains(">readA some desc"));
    assert!(out.contains(">readB"));
}

#[test]
fn bad_orientation_errors() {
    let td = tempfile::tempdir().unwrap();
    let fasta_p = td.path().join("in.fa");
    let tsv_p = td.path().join("map.tsv");
    write(&fasta_p, FASTA);
    write(&tsv_p, BAD_TSV);

    let mut cmd = Command::cargo_bin("restrand-fasta").unwrap();
    cmd.args([
        "-f",
        fasta_p.to_str().unwrap(),
        "-t",
        tsv_p.to_str().unwrap(),
    ]);
    cmd.assert()
        .failure()
        .stderr(predicate::str::contains("Unrecognized orientation value"));
}
