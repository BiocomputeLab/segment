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
const SEG_B: &str = "TTGCAACCGTAAGGTCTTGA";
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
        .args(["--classifications", output.to_str().unwrap()])
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
        .args(["--classifications", output.to_str().unwrap()])
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

/// --per-segment-scores reads each segment's threshold from square brackets after its
/// name, and the brackets do not reach the output. The same file without the flag keeps
/// the brackets as part of the name and holds every segment to --min-norm-score, so a
/// file written before the option existed is unaffected by it.
#[test]
fn cli_per_segment_scores_reads_thresholds_from_segment_names() {
    let dir = TempDir::new().unwrap();
    let refs = dir.path().join("refs.fasta");
    // A is planted with one mismatch, so it scores 37/20 = 1.85: above 1.8, below 1.9.
    let degraded_a = "CGATGCTAGCAACGGATCAT";
    assert_eq!(degraded_a.len(), SEG_A.len());
    let sequences = dir.path().join("reads.fastq");
    write_fastq(
        &sequences,
        &[("read", &format!("{JUNK}{START}{degraded_a}{END}{JUNK}"))],
    );

    let classify = |threshold: &str, extra: &[&str]| -> String {
        fs::write(
            &refs,
            format!(">start\n{START}\n>end\n{END}\n>A[{threshold}]\n{SEG_A}\n"),
        )
        .unwrap();
        let out = dir.path().join("out.txt");
        let mut args = vec!["--start-end-segs"];
        args.extend_from_slice(extra);
        run_segment(&refs, &sequences, &out, &args);
        fs::read_to_string(&out).unwrap()
    };

    // With the flag, the bracketed score decides whether A is reported, and the name it
    // is reported under has no brackets in it.
    assert_eq!(
        classify("1.8", &["--per-segment-scores"]),
        "read\tstart-A-end\n"
    );
    assert_eq!(
        classify("1.9", &["--per-segment-scores"]),
        "read\tstart-end\n"
    );

    // Without it the name is still reported without its brackets, but the score they
    // hold is ignored and the global threshold of 1.5 applies instead - which the
    // degraded copy clears, so A is reported after all.
    assert_eq!(classify("1.9", &[]), "read\tstart-A-end\n");
}

/// Writing a score and forgetting the flag silently gets you the global threshold
/// instead of the one you asked for, so the run says so rather than leaving it to be
/// discovered from the results.
#[test]
fn cli_warns_when_a_segment_names_a_score_without_the_flag() {
    let dir = TempDir::new().unwrap();
    let refs = dir.path().join("refs.fasta");
    fs::write(&refs, format!(">A[1.9]\n{SEG_A}\n>B\n{SEG_B}\n")).unwrap();

    let sequences = dir.path().join("reads.fastq");
    write_fastq(&sequences, &[("read", &format!("{JUNK}{SEG_A}{JUNK}"))]);

    let output = dir.path().join("out.txt");
    let result = Command::new(env!("CARGO_BIN_EXE_segment"))
        .args(["--segments", refs.to_str().unwrap()])
        .args(["--sequences", sequences.to_str().unwrap()])
        .args(["--classifications", output.to_str().unwrap()])
        .args(["--min-seq-len", "0", "--max-seq-len", "100000"])
        .output()
        .unwrap();
    assert!(result.status.success(), "it is a warning, not an error");
    let stderr = String::from_utf8(result.stderr).unwrap();
    assert!(
        stderr.contains("segment 'A[1.9]'") && stderr.contains("--per-segment-scores"),
        "{stderr}"
    );
    // The segment named no score, so nothing is said about it.
    assert!(!stderr.contains("segment 'B'"), "{stderr}");
    // ...and the results use the stripped name regardless.
    assert_eq!(fs::read_to_string(&output).unwrap(), "read\tA\n");
}

/// A bracketed value that is not a usable score stops the run and names the segment,
/// rather than leaving a segment that is silently never found.
#[test]
fn cli_rejects_an_unreadable_per_segment_score() {
    let dir = TempDir::new().unwrap();
    let refs = dir.path().join("refs.fasta");
    fs::write(&refs, format!(">A[high]\n{SEG_A}\n")).unwrap();

    let sequences = dir.path().join("reads.fastq");
    write_fastq(&sequences, &[("read", &format!("{JUNK}{SEG_A}{JUNK}"))]);

    let result = Command::new(env!("CARGO_BIN_EXE_segment"))
        .args(["--segments", refs.to_str().unwrap()])
        .args(["--sequences", sequences.to_str().unwrap()])
        .args([
            "--classifications",
            dir.path().join("out.txt").to_str().unwrap(),
        ])
        .arg("--per-segment-scores")
        .output()
        .unwrap();
    assert!(!result.status.success(), "the run should have failed");
    let stderr = String::from_utf8(result.stderr).unwrap();
    assert!(stderr.contains("Segment 'A[high]'"), "{stderr}");
    assert!(stderr.contains("(record 1)"), "{stderr}");
    assert!(stderr.contains("'high' is not a number"), "{stderr}");
    assert!(
        stderr.contains("'SEG1[1.5]'"),
        "the syntax is shown: {stderr}"
    );
}

/// --counts writes a CSV of every distinct classification and how many reads produced
/// it, alongside the per-read results rather than instead of them. The counts must
/// reconcile with the results file exactly, since a tally that does not add up to the
/// reads is worse than no tally at all.
#[test]
fn cli_counts_flag_writes_classification_counts() {
    let dir = TempDir::new().unwrap();
    let refs = dir.path().join("refs.fasta");
    fs::write(
        &refs,
        format!(">start\n{START}\n>end\n{END}\n>A\n{SEG_A}\n"),
    )
    .unwrap();

    let with_a = format!("{JUNK}{START}{SEG_A}{END}{JUNK}");
    let without_a = format!("{JUNK}{START}{END}{JUNK}");
    let sequences = dir.path().join("reads.fastq");
    write_fastq(
        &sequences,
        &[
            ("a1", &with_a),
            ("bare", &without_a),
            ("a2", &rc(&with_a)),
            ("dropped", &format!("{JUNK}{SEG_A}{JUNK}")),
            ("a3", &with_a),
        ],
    );

    let output = dir.path().join("out.txt");
    let counts = dir.path().join("counts.csv");
    run_segment(
        &refs,
        &sequences,
        &output,
        &["--start-end-segs", "--counts", counts.to_str().unwrap()],
    );

    // Most frequent first. The read with no anchors never reaches a classification, so
    // it is absent here and reported in the run summary on stderr instead.
    assert_eq!(
        fs::read_to_string(&counts).unwrap(),
        "segments,count\nstart-A-end,3\nstart-end,1\n"
    );
    // The per-read results are still written, and the counts add up to them.
    let results = fs::read_to_string(&output).unwrap();
    assert_eq!(results.lines().count(), 4);
    assert_eq!(
        results
            .lines()
            .filter(|l| l.ends_with("start-A-end"))
            .count(),
        3
    );
}

/// Without --counts no CSV is written at all, so the flag adds a file rather than
/// changing what a run without it produces.
#[test]
fn cli_writes_no_counts_csv_unless_asked() {
    let dir = TempDir::new().unwrap();
    let refs = dir.path().join("refs.fasta");
    fs::write(
        &refs,
        format!(">start\n{START}\n>end\n{END}\n>A\n{SEG_A}\n"),
    )
    .unwrap();

    let sequences = dir.path().join("reads.fastq");
    write_fastq(
        &sequences,
        &[("read", &format!("{JUNK}{START}{SEG_A}{END}{JUNK}"))],
    );

    let plain = dir.path().join("plain.txt");
    run_segment(&refs, &sequences, &plain, &["--start-end-segs"]);
    assert_eq!(fs::read_to_string(&plain).unwrap(), "read\tstart-A-end\n");
    assert!(!dir.path().join("counts.csv").exists());

    // ...and asking for one does not change the per-read results.
    let with_counts = dir.path().join("with_counts.txt");
    let counts = dir.path().join("counts.csv");
    run_segment(
        &refs,
        &sequences,
        &with_counts,
        &["--start-end-segs", "--counts", counts.to_str().unwrap()],
    );
    assert_eq!(
        fs::read_to_string(&with_counts).unwrap(),
        fs::read_to_string(&plain).unwrap()
    );
    assert_eq!(
        fs::read_to_string(&counts).unwrap(),
        "segments,count\nstart-A-end,1\n"
    );
}

/// An unwritable --counts path is reported before the reads are processed rather than
/// after, so a long run does not finish only to find it cannot write the file.
#[test]
fn cli_rejects_an_unwritable_counts_path_up_front() {
    let dir = TempDir::new().unwrap();
    let refs = dir.path().join("refs.fasta");
    fs::write(
        &refs,
        format!(">start\n{START}\n>end\n{END}\n>A\n{SEG_A}\n"),
    )
    .unwrap();

    let sequences = dir.path().join("reads.fastq");
    write_fastq(
        &sequences,
        &[("read", &format!("{JUNK}{START}{SEG_A}{END}{JUNK}"))],
    );

    let output = dir.path().join("out.txt");
    let result = Command::new(env!("CARGO_BIN_EXE_segment"))
        .args(["--segments", refs.to_str().unwrap()])
        .args(["--sequences", sequences.to_str().unwrap()])
        .args(["--classifications", output.to_str().unwrap()])
        // A directory that does not exist, so the file cannot be created.
        .args([
            "--counts",
            dir.path().join("nope/counts.csv").to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(!result.status.success(), "the run should have failed");
    let stderr = String::from_utf8(result.stderr).unwrap();
    assert!(stderr.contains("Could not create counts file"), "{stderr}");
}

/// A degenerate segment reaches the classifier through the whole pipeline: one segment
/// carrying an R (A or G) finds both reads it stands for, on either strand, and the
/// warning about ambiguity is not raised for a segment this specific.
#[test]
fn cli_finds_both_sequences_a_degenerate_segment_stands_for() {
    let dir = TempDir::new().unwrap();
    let refs = dir.path().join("refs.fasta");
    // SEG_A with an R where it has an A, so it also matches the G-carrying variant.
    let degenerate = "CGRTGCTAGCTACGGATCAT";
    let variant = "CGGTGCTAGCTACGGATCAT";
    fs::write(
        &refs,
        format!(">start\n{START}\n>end\n{END}\n>A\n{degenerate}\n"),
    )
    .unwrap();

    let sequences = dir.path().join("reads.fastq");
    let spelled_out = format!("{JUNK}{START}{SEG_A}{END}{JUNK}");
    let with_variant = format!("{JUNK}{START}{variant}{END}{JUNK}");
    write_fastq(
        &sequences,
        &[
            ("exact", &spelled_out),
            ("variant", &with_variant),
            ("variant_rc", &rc(&with_variant)),
        ],
    );

    let output = dir.path().join("out.txt");
    let result = Command::new(env!("CARGO_BIN_EXE_segment"))
        .args(["--segments", refs.to_str().unwrap()])
        .args(["--sequences", sequences.to_str().unwrap()])
        .args(["--classifications", output.to_str().unwrap()])
        .args(["--min-seq-len", "0", "--max-seq-len", "100000"])
        .args(["--min-norm-score", "1.5", "--start-end-segs"])
        .output()
        .unwrap();
    assert!(result.status.success());
    assert_eq!(
        fs::read_to_string(&output).unwrap(),
        "exact\tstart-A-end\nvariant\tstart-A-end\nvariant_rc\tstart-A-end\n"
    );
    let stderr = String::from_utf8(result.stderr).unwrap();
    assert!(
        !stderr.contains("ambiguous enough"),
        "one degenerate base should not trigger the ambiguity warning:\n{stderr}"
    );
}

/// A segment degenerate enough to clear the threshold on random sequence is not an
/// error - it is a legitimate, if unusual, thing to search for - but it must be called
/// out by name so nobody mistakes the resulting flood of hits for a real result.
#[test]
fn cli_warns_about_a_segment_that_would_match_anything() {
    let dir = TempDir::new().unwrap();
    let refs = dir.path().join("refs.fasta");
    // Sixteen Ns in twenty bases matches 85% of random sequence, above the 83% the
    // default threshold demands.
    let vague = format!("{}{}", "N".repeat(16), &SEG_A[16..]);
    fs::write(&refs, format!(">A\n{SEG_A}\n>vague\n{vague}\n")).unwrap();

    let sequences = dir.path().join("reads.fastq");
    write_fastq(&sequences, &[("read", &format!("{JUNK}{SEG_A}{JUNK}"))]);

    let output = dir.path().join("out.txt");
    let result = Command::new(env!("CARGO_BIN_EXE_segment"))
        .args(["--segments", refs.to_str().unwrap()])
        .args(["--sequences", sequences.to_str().unwrap()])
        .args(["--classifications", output.to_str().unwrap()])
        .args(["--min-seq-len", "0", "--max-seq-len", "100000"])
        .args(["--min-norm-score", "1.5"])
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "a vague segment is a warning, not an error"
    );
    let stderr = String::from_utf8(result.stderr).unwrap();
    assert!(
        stderr.contains("segment 'vague' is ambiguous enough"),
        "the vague segment should be named:\n{stderr}"
    );
    assert!(
        !stderr.contains("segment 'A'"),
        "the specific segment should not be named:\n{stderr}"
    );
}

/// A segment carrying a character that is no nucleotide could never be found, so the
/// run stops with a message naming the segment and the offending character rather than
/// quietly reporting that it is absent from every read.
#[test]
fn cli_rejects_a_segment_that_is_not_a_nucleotide_sequence() {
    let dir = TempDir::new().unwrap();
    let refs = dir.path().join("refs.fasta");
    fs::write(&refs, format!(">A\n{SEG_A}\n>bad\nACGTXACGT\n")).unwrap();

    let sequences = dir.path().join("reads.fastq");
    write_fastq(&sequences, &[("read", &format!("{JUNK}{SEG_A}{JUNK}"))]);

    let output = dir.path().join("out.txt");
    let result = Command::new(env!("CARGO_BIN_EXE_segment"))
        .args(["--segments", refs.to_str().unwrap()])
        .args(["--sequences", sequences.to_str().unwrap()])
        .args(["--classifications", output.to_str().unwrap()])
        .output()
        .unwrap();
    assert!(!result.status.success(), "the run should have failed");
    let stderr = String::from_utf8(result.stderr).unwrap();
    assert!(stderr.contains("Segment 'bad'"), "{stderr}");
    assert!(stderr.contains("'X' at position 5"), "{stderr}");
    assert!(!output.exists(), "no output should have been written");
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
