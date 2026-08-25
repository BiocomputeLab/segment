// End-to-end tests driving the compiled `segment` binary. The unit tests in
// src/tests.rs cover the classification logic; these check the user-facing contract:
// argument handling, format auto-detection on a real file, and the output written to
// disk.

use flate2::Compression;
use flate2::write::GzEncoder;
use noodles_bam as bam;
use noodles_sam::{self as sam, alignment::io::Write as _, alignment::record::Flags};
use std::fs::{self, File};
use std::io::Write as _;
use std::path::Path;
use std::process::Command;
use tempfile::TempDir;

const START: &str = "AGGCATTCGAGCTTAACGGT";
const END: &str = "TCCAGTTAGCCTGAAGTCAC";
const SEG_A: &str = "CGATGCTAGCTACGGATCAT";
const JUNK: &str = "GGTTAACCGG";

/// Reverse complement of a sequence.
fn rc(seq: &str) -> String {
    seq.chars()
        .rev()
        .map(|c| match c {
            'A' => 'T',
            'T' => 'A',
            'C' => 'G',
            'G' => 'C',
            other => other,
        })
        .collect()
}

/// Write a FASTQ, giving every base the same placeholder quality.
fn write_fastq(path: &Path, reads: &[(&str, &str)]) {
    let mut f = File::create(path).unwrap();
    for (name, seq) in reads {
        writeln!(f, "@{}\n{}\n+\n{}", name, seq, "I".repeat(seq.len())).unwrap();
    }
}

/// Write a gzipped FASTQ.
fn write_fastq_gz(path: &Path, reads: &[(&str, &str)]) {
    let mut encoder = GzEncoder::new(File::create(path).unwrap(), Compression::default());
    for (name, seq) in reads {
        writeln!(encoder, "@{}\n{}\n+\n{}", name, seq, "I".repeat(seq.len())).unwrap();
    }
    encoder.finish().unwrap();
}

/// Write an unaligned BAM.
fn write_bam(path: &Path, reads: &[(&str, &str)]) {
    let mut writer = bam::io::Writer::new(File::create(path).unwrap());
    let header = sam::Header::default();
    writer.write_alignment_header(&header).unwrap();
    for (name, seq) in reads {
        let record = sam::alignment::RecordBuf::builder()
            .set_name(*name)
            .set_flags(Flags::UNMAPPED)
            .set_sequence(seq.as_bytes().to_vec().into())
            .build();
        writer.write_alignment_record(&header, &record).unwrap();
    }
    writer.finish(&header).unwrap();
}

/// Run the compiled binary with fixed bounds and threshold, asserting it succeeds.
fn run_segment(refs: &Path, sequences: &Path, output: &Path, extra: &[&str]) {
    let status = Command::new(env!("CARGO_BIN_EXE_segment"))
        .args(["--segments", refs.to_str().unwrap()])
        .args(["--sequences", sequences.to_str().unwrap()])
        .args(["--output", output.to_str().unwrap()])
        .args(["--min-seq-len", "0", "--max-seq-len", "100000"])
        .args(["--min-norm-score", "1.5"])
        .args(extra)
        .status()
        .unwrap();
    assert!(status.success(), "segment exited with {status}");
}

/// The CLI accepts FASTQ, gzipped FASTQ and unaligned BAM interchangeably with no
/// format flag and no filename convention, and writes the identical tab-separated
/// "<read name>\t<segments>" output from each. The two reads describe the same
/// construct on opposite strands, so both must resolve to the same segment string.
#[test]
fn cli_accepts_every_input_format_and_writes_matching_output() {
    let dir = TempDir::new().unwrap();
    let refs = dir.path().join("refs.fasta");
    fs::write(
        &refs,
        format!(">start\n{START}\n>end\n{END}\n>A\n{SEG_A}\n"),
    )
    .unwrap();

    let forward = format!("{JUNK}{START}{SEG_A}{END}{JUNK}");
    let reverse = rc(&forward);
    let reads = [("fwd", forward.as_str()), ("rev", reverse.as_str())];

    let fastq = dir.path().join("reads.fastq");
    let gzipped = dir.path().join("reads.fastq.gz");
    let bam_file = dir.path().join("reads.bam");
    write_fastq(&fastq, &reads);
    write_fastq_gz(&gzipped, &reads);
    write_bam(&bam_file, &reads);

    for input in [&fastq, &gzipped, &bam_file] {
        let output = dir.path().join("out.txt");
        run_segment(&refs, input, &output, &["--start-end-segs"]);
        assert_eq!(
            fs::read_to_string(&output).unwrap(),
            "fwd\tstart-A-end\nrev\tstart-A-end\n",
            "unexpected output for input {input:?}"
        );
    }
}

/// The documented default score threshold is what the binary actually applies. Pinning
/// it here means the README and the CLI cannot drift apart unnoticed.
#[test]
fn cli_default_score_threshold_is_documented_value() {
    let out = Command::new(env!("CARGO_BIN_EXE_segment"))
        .arg("--help")
        .output()
        .unwrap();
    let help = String::from_utf8(out.stdout).unwrap();
    assert!(
        help.contains("[default: 1.5]"),
        "--min-norm-score default is not 1.5:\n{help}"
    );
}

/// The summary is written to stderr, never to stdout or the results file, so it cannot
/// contaminate a downstream pipeline reading either.
#[test]
fn cli_writes_summary_to_stderr_only() {
    let dir = TempDir::new().unwrap();
    let refs = dir.path().join("refs.fasta");
    fs::write(&refs, format!(">start\n{START}\n>end\n{END}\n>A\n{SEG_A}\n")).unwrap();

    let sequences = dir.path().join("reads.fastq");
    write_fastq(
        &sequences,
        &[
            ("good", &format!("{JUNK}{START}{SEG_A}{END}{JUNK}")),
            ("no_anchors", &format!("{JUNK}{SEG_A}{JUNK}")),
        ],
    );

    let output = dir.path().join("out.txt");
    let result = Command::new(env!("CARGO_BIN_EXE_segment"))
        .args(["--segments", refs.to_str().unwrap()])
        .args(["--sequences", sequences.to_str().unwrap()])
        .args(["--output", output.to_str().unwrap()])
        .args(["--min-seq-len", "0", "--max-seq-len", "100000"])
        .args(["--min-norm-score", "1.5", "--start-end-segs"])
        .output()
        .unwrap();
    assert!(result.status.success());

    let stderr = String::from_utf8(result.stderr).unwrap();
    assert!(stderr.contains("Summary"), "no summary on stderr:\n{stderr}");
    assert!(stderr.contains("classified:"), "{stderr}");
    assert!(
        stderr.contains("neither segment found"),
        "the unanchored read should be reported as a rejection:\n{stderr}"
    );

    assert!(String::from_utf8(result.stdout).unwrap().is_empty(), "stdout must stay clean");
    // The results file holds only the classified read, with no summary text in it.
    assert_eq!(fs::read_to_string(&output).unwrap(), "good\tstart-A-end\n");
}

/// The --circular flag reaches the pipeline: a read that wraps the origin (END before
/// START) is dropped without it and recovered with it, from the same input file.
#[test]
fn cli_circular_flag_recovers_reads_that_wrap_the_origin() {
    let dir = TempDir::new().unwrap();
    let refs = dir.path().join("refs.fasta");
    fs::write(
        &refs,
        format!(">start\n{START}\n>end\n{END}\n>A\n{SEG_A}\n"),
    )
    .unwrap();

    let sequences = dir.path().join("reads.fastq");
    write_fastq(&sequences, &[("wrapped", &format!("{END}{START}{SEG_A}"))]);

    let linear = dir.path().join("linear.txt");
    run_segment(&refs, &sequences, &linear, &["--start-end-segs"]);
    assert_eq!(fs::read_to_string(&linear).unwrap(), "");

    let circular = dir.path().join("circular.txt");
    run_segment(&refs, &sequences, &circular, &["--start-end-segs", "--circular"]);
    assert_eq!(
        fs::read_to_string(&circular).unwrap(),
        "wrapped\tstart-A-end\n"
    );
}
