// Test suite for the segment tool.
//
// Every fixture is deliberately small so that the expected result can be derived by
// hand. All segments are 20 bp and mutually dissimilar, and all tests use a minimum
// normalised score of 1.5, which for an ungapped 20 bp segment demands >= 17/20
// matching bases. A segment therefore only aligns above threshold at the position it
// was deliberately planted, which keeps every expectation below unambiguous.

use super::*;

use flate2::Compression;
use flate2::write::GzEncoder;
use noodles_sam::{self as sam, alignment::io::Write as _, alignment::record::Flags};
use tempfile::TempDir;

// ---------------------------------------------------------------------------
// Fixtures and helpers
// ---------------------------------------------------------------------------

const START: &str = "AGGCATTCGAGCTTAACGGT";
const END: &str = "TCCAGTTAGCCTGAAGTCAC";
const SEG_A: &str = "CGATGCTAGCTACGGATCAT";
const SEG_B: &str = "TTGCAACCGTAAGGTCTTGA";
/// Filler placed between segments; matches nothing above threshold.
const SPACER: &str = "TTAGGACCTGAT";
/// A segment whose first six bases repeat as its last six ("CGATGC"), so two copies
/// placed six bases apart agree on the bases they share. This is what makes it possible
/// for a repeated segment to overlap itself by more than the 3 bp tolerance while both
/// copies still score above the threshold.
const SEG_REPEAT: &str = "CGATGCTAGCTACGCGATGC";
/// Adapter-like junk outside the start/end anchors.
const JUNK_5: &str = "GGTTAACCGG";
const JUNK_3: &str = "CCAATTGGCC";

const MIN_SCORE: f32 = 1.5;

/// Reverse complement of a sequence.
fn rc(seq: &str) -> String {
    String::from_utf8(dna::revcomp(seq.as_bytes())).unwrap()
}

/// Transverse the base at each given offset (A<->C, G<->T), introducing exactly
/// `positions.len()` mismatches at known places.
fn mutate(seq: &str, positions: &[usize]) -> String {
    let mut bases = seq.as_bytes().to_vec();
    for &p in positions {
        bases[p] = match bases[p] {
            b'A' => b'C',
            b'C' => b'A',
            b'G' => b'T',
            _ => b'G',
        };
    }
    String::from_utf8(bases).unwrap()
}

/// Build the name-to-sequence map the classifier takes.
fn seg_map(pairs: &[(&str, &str)]) -> HashMap<String, Segment> {
    seg_map_scored(
        &pairs
            .iter()
            .map(|&(n, s)| (n, s, MIN_SCORE))
            .collect::<Vec<_>>(),
    )
}

/// The same, with an explicit minimum normalised score per segment, for the tests that
/// care which threshold each one is held to.
fn seg_map_scored(triples: &[(&str, &str, f32)]) -> HashMap<String, Segment> {
    triples
        .iter()
        .map(|(name, seq, min_norm_score)| {
            (
                name.to_string(),
                Segment {
                    seq: seq.as_bytes().to_vec(),
                    min_norm_score: *min_norm_score,
                },
            )
        })
        .collect()
}

/// Load a segments file at the default threshold, with per-segment scores switched off.
fn load(path: &Path) -> Result<HashMap<String, Segment>> {
    load_fasta(path, MIN_SCORE, false)
}

/// ...and with them switched on.
fn load_scored(path: &Path) -> Result<HashMap<String, Segment>> {
    load_fasta(path, MIN_SCORE, true)
}

/// Classify a read exactly as given, with no start/end anchors supplied.
fn classify(segments: &[(&str, &str)], read: &str) -> String {
    segment_string(&classify_read_segments(
        &seg_map(segments),
        read.as_bytes(),
        false,
    ))
}

/// Classify a read that has already been trimmed to run from 'start' to 'end'.
fn classify_anchored(segments: &[(&str, &str)], read: &str) -> String {
    // Zero-length anchor spans: these fixtures are the interior of a trimmed read, so
    // the anchors bound it without occupying any of it. `classified_read` is what puts
    // them back, so going through it keeps the helper honest about where they come from.
    classified_read(
        "read".to_string(),
        read.as_bytes(),
        &seg_map(segments),
        Some((0, 0)),
        false,
    )
    .segments
}

/// `cut_reorient_seq` over strings, for readable expectations.
fn cut(seq: &str, start: usize, end: usize, revcomp: bool) -> String {
    String::from_utf8(cut_reorient_seq(seq.as_bytes(), start, end, revcomp)).unwrap()
}

/// Align one segment against a read with the same scoring the classifier uses
/// (+2 match, -1 mismatch, -2 gap).
fn align(segment: &str, read: &str, min_score: i32) -> Vec<MultiAlignment> {
    MultiTracebackAligner::new(2, -1, -2).align_multiple(
        "seg",
        segment.as_bytes(),
        read.as_bytes(),
        min_score,
    )
}

/// Unwrap the single expected alignment as (score, normalised score, start, end).
fn only(alignments: Vec<MultiAlignment>) -> (i32, f32, usize, usize) {
    assert_eq!(
        alignments.len(),
        1,
        "expected exactly one alignment, got {alignments:?}"
    );
    let a = &alignments[0];
    (a.score, a.norm_score, a.target_start, a.target_end)
}

/// A read as the pipeline expects it.
fn raw(name: &str, seq: &str) -> RawRead {
    RawRead {
        name: name.to_string(),
        seq: seq.as_bytes().to_vec(),
    }
}

/// Reads as comparable (name, sequence) pairs.
fn as_pairs(reads: &[RawRead]) -> Vec<(String, String)> {
    reads
        .iter()
        .map(|r| (r.name.clone(), String::from_utf8(r.seq.clone()).unwrap()))
        .collect()
}

/// Run the full per-read pipeline with read-length filtering disabled.
fn run(
    reads: Vec<RawRead>,
    segments: &[(&str, &str)],
    circular: bool,
    start_end_segs: bool,
) -> Vec<(String, String)> {
    process_reads(
        reads,
        &seg_map(segments),
        (0, 0),
        circular,
        start_end_segs,
        false,
    )
    .unwrap()
    .0
    .into_iter()
    .map(|r| (r.name, r.segments))
    .collect()
}

/// The same with --detailed-output on, as (segments, located, extracted sequence) for
/// the single read the caller passed in.
fn run_detailed(
    read: &str,
    segments: &[(&str, &str)],
    circular: bool,
    start_end_segs: bool,
) -> (String, String, String) {
    let (mut classified, _) = process_reads(
        vec![raw("read", read)],
        &seg_map(segments),
        (0, 0),
        circular,
        start_end_segs,
        true,
    )
    .unwrap();
    assert_eq!(classified.len(), 1, "the read should have been classified");
    let r = classified.remove(0);
    let detail = r
        .detail
        .expect("--detailed-output should have filled the detail");
    (r.segments, detail.located, detail.sequence)
}

/// The full construct in sequencing orientation: START-A-SPACER-B-END.
fn construct() -> String {
    format!("{}{}{}{}{}", START, SEG_A, SPACER, SEG_B, END)
}

/// The segment set for anchored runs: both anchors plus A and B.
fn anchored_segments() -> Vec<(&'static str, &'static str)> {
    vec![("start", START), ("end", END), ("A", SEG_A), ("B", SEG_B)]
}

/// The expected result for a single classified read.
fn expect(name: &str, segments: &str) -> Vec<(String, String)> {
    vec![(name.to_string(), segments.to_string())]
}

/// Assert a call failed, and that its message mentions `needle`.
fn assert_error<T>(result: Result<T>, needle: &str) {
    match result {
        Ok(_) => panic!("expected an error mentioning {needle:?}, but the call succeeded"),
        Err(e) => assert!(
            e.contains(needle),
            "error {e:?} does not mention {needle:?}"
        ),
    }
}

/// Write a segments FASTA.
fn write_fasta(path: &Path, records: &[(&str, &str)]) {
    let mut f = File::create(path).unwrap();
    for (name, seq) in records {
        writeln!(f, ">{}\n{}", name, seq).unwrap();
    }
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

/// Write an unaligned BAM. `seq` is stored verbatim, so setting
/// `Flags::REVERSE_COMPLEMENTED` declares it to be the reverse complement of the read.
fn write_bam(path: &Path, reads: &[(&str, &str, Flags)]) {
    let mut writer = bam::io::Writer::new(File::create(path).unwrap());
    let header = sam::Header::default();
    writer.write_alignment_header(&header).unwrap();
    for (name, seq, flags) in reads {
        let record = sam::alignment::RecordBuf::builder()
            .set_name(*name)
            .set_flags(*flags)
            .set_sequence(seq.as_bytes().to_vec().into())
            .build();
        writer.write_alignment_record(&header, &record).unwrap();
    }
    writer.finish(&header).unwrap();
}

// ---------------------------------------------------------------------------
// Input format detection and loading
// ---------------------------------------------------------------------------

/// The same reads stored as plain FASTQ, gzipped FASTQ and unaligned BAM must each be
/// detected from their contents alone (no extension hints, no flags) and yield
/// identical names and sequences, so results cannot depend on the container used.
#[test]
fn all_input_formats_are_detected_and_load_identical_reads() {
    let dir = TempDir::new().unwrap();
    let reads = [("r1", SEG_A), ("r2", SEG_B), ("r3", SPACER)];

    let fastq = dir.path().join("reads.fastq");
    let gzipped = dir.path().join("reads.fastq.gz");
    let bam_file = dir.path().join("reads.bam");
    write_fastq(&fastq, &reads);
    write_fastq_gz(&gzipped, &reads);
    let bam_reads: Vec<_> = reads
        .iter()
        .map(|(n, s)| (*n, *s, Flags::UNMAPPED))
        .collect();
    write_bam(&bam_file, &bam_reads);

    assert!(matches!(detect_format(&fastq).unwrap(), InputFormat::Fastq));
    assert!(matches!(
        detect_format(&gzipped).unwrap(),
        InputFormat::FastqGz
    ));
    // BAM is itself gzip (BGZF) framed, so detection must look past the shared gzip
    // magic bytes at the decompressed "BAM\1" magic to tell it from a gzipped FASTQ.
    assert!(matches!(
        detect_format(&bam_file).unwrap(),
        InputFormat::Bam
    ));

    let expected: Vec<_> = reads
        .iter()
        .map(|(n, s)| (n.to_string(), s.to_string()))
        .collect();
    assert_eq!(as_pairs(&load_reads(&fastq).unwrap().reads), expected);
    assert_eq!(as_pairs(&load_reads(&gzipped).unwrap().reads), expected);
    assert_eq!(as_pairs(&load_reads(&bam_file).unwrap().reads), expected);
}

/// BAM stores the sequence of a reverse-strand record already reverse complemented.
/// Loading must undo that so every read reaches the classifier in the orientation it
/// came off the sequencer, matching what the equivalent FASTQ would have contained.
#[test]
fn bam_reverse_complemented_records_are_restored_to_read_orientation() {
    let dir = TempDir::new().unwrap();
    let bam_file = dir.path().join("reads.bam");
    let stored = rc(SEG_A);
    write_bam(
        &bam_file,
        &[
            ("fwd", SEG_A, Flags::UNMAPPED),
            (
                "rev",
                stored.as_str(),
                Flags::UNMAPPED | Flags::REVERSE_COMPLEMENTED,
            ),
        ],
    );

    assert_eq!(
        as_pairs(&load_reads(&bam_file).unwrap().reads),
        vec![
            ("fwd".to_string(), SEG_A.to_string()),
            ("rev".to_string(), SEG_A.to_string()),
        ]
    );
}

/// A file that is neither FASTQ nor gzip framed is refused with a message naming the
/// file and the formats accepted - an error value, not a panic.
#[test]
fn unrecognised_input_is_rejected() {
    let dir = TempDir::new().unwrap();
    let path = dir.path().join("reads.txt");
    File::create(&path)
        .unwrap()
        .write_all(b">this is a fasta\nACGT\n")
        .unwrap();
    assert_error(detect_format(&path), "Could not determine the format of");
    assert_error(detect_format(&path), "reads.txt");
}

/// Gzip-framed input that decompresses to neither FASTQ nor BAM is refused as an
/// unusable compressed file, distinctly from the uncompressed case above.
#[test]
fn unrecognised_gzipped_input_is_rejected() {
    let dir = TempDir::new().unwrap();
    let path = dir.path().join("reads.gz");
    let mut encoder = GzEncoder::new(File::create(&path).unwrap(), Compression::default());
    encoder.write_all(b">this is a fasta\nACGT\n").unwrap();
    encoder.finish().unwrap();
    assert_error(detect_format(&path), "gzip-compressed");
}

/// A missing input file names the file rather than panicking somewhere inside the
/// FASTA/FASTQ machinery.
#[test]
fn missing_files_are_reported_by_name() {
    let dir = TempDir::new().unwrap();
    let missing = dir.path().join("nope.fasta");
    assert_error(load(&missing), "Could not open segments file");
    assert_error(detect_format(&missing), "Could not open sequences file");
}

/// Segment definitions are keyed by FASTA record ID (the first word of the header),
/// which is what later stages use both to find "start"/"end" and to label results.
#[test]
fn segment_definitions_are_loaded_from_fasta() {
    let dir = TempDir::new().unwrap();
    let path = dir.path().join("refs.fasta");
    let mut f = File::create(&path).unwrap();
    writeln!(f, ">start some description\n{}\n>A\n{}", START, SEG_A).unwrap();
    drop(f);

    let segments = load(&path).unwrap();
    assert_eq!(segments.len(), 2);
    assert_eq!(segments["start"].seq, START.as_bytes());
    assert_eq!(segments["A"].seq, SEG_A.as_bytes());
}

// ---------------------------------------------------------------------------
// Input validation and recovery
// ---------------------------------------------------------------------------

/// A segment record with no sequence used to divide by its own length, producing a NaN
/// score that crashed the sort with a bare "called `Option::unwrap()` on a `None`
/// value". It is now refused up front, naming the offending segment.
#[test]
fn empty_segment_sequence_is_rejected() {
    let dir = TempDir::new().unwrap();
    let path = dir.path().join("refs.fasta");
    write_fasta(&path, &[("start", START), ("broken", ""), ("A", SEG_A)]);
    assert_error(load(&path), "Segment 'broken'");
    assert_error(load(&path), "has no sequence");
}

/// Defence in depth for the same bug: if a zero-length segment reached the classifier
/// some other way it must be ignored - not crash the score sort with a NaN, and not
/// "match" at every position on a score of zero.
#[test]
fn zero_length_segment_does_not_crash_the_classifier() {
    let read = format!("{}{}", SEG_A, SPACER);
    assert_eq!(classify(&[("A", SEG_A), ("empty", "")], &read), "A");
}

/// A character outside the nucleotide alphabet matches nothing, so a segment carrying
/// one could never be found. Left alone that is a silently empty result; it is refused
/// instead, naming the character and where it is.
#[test]
fn non_nucleotide_segment_characters_are_rejected() {
    let dir = TempDir::new().unwrap();
    let path = dir.path().join("refs.fasta");
    write_fasta(&path, &[("A", SEG_A), ("bad", "ACGTXACGT")]);
    assert_error(load(&path), "Segment 'bad'");
    assert_error(load(&path), "has 'X' at position 5");
}

/// A segment of nothing but N matches perfectly at every position of every read, which
/// is never a question anyone meant to ask, so it is refused up front.
#[test]
fn entirely_ambiguous_segment_is_rejected() {
    let dir = TempDir::new().unwrap();
    let path = dir.path().join("refs.fasta");
    write_fasta(&path, &[("A", SEG_A), ("blank", "NNNNNNNNNN")]);
    assert_error(load(&path), "Segment 'blank'");
    assert_error(load(&path), "entirely ambiguous");
}

/// A partly ambiguous segment is fine, though - only every base being ambiguous is not.
#[test]
fn partly_ambiguous_segment_is_accepted() {
    let dir = TempDir::new().unwrap();
    let path = dir.path().join("refs.fasta");
    write_fasta(&path, &[("A", "NNNNNCGATGCNNNNN")]);
    assert_eq!(load(&path).unwrap()["A"].seq, b"NNNNNCGATGCNNNNN");
}

/// An RNA spelling of a segment should find the same reads as the DNA one. U is folded
/// into T on load rather than in the aligner, because U has no complement of its own and
/// would otherwise survive reverse complementing unchanged.
#[test]
fn uracil_in_a_segment_is_read_as_thymine() {
    let dir = TempDir::new().unwrap();
    let path = dir.path().join("refs.fasta");
    let rna = SEG_A.replace('T', "U");
    write_fasta(&path, &[("A", &rna)]);
    assert_eq!(load(&path).unwrap()["A"].seq, SEG_A.as_bytes());
}

/// Two segments sharing a name make it ambiguous which sequence to search for, so
/// processing stops rather than silently keeping whichever was read last.
#[test]
fn duplicate_segment_names_are_rejected() {
    let dir = TempDir::new().unwrap();
    let path = dir.path().join("refs.fasta");
    write_fasta(&path, &[("A", SEG_A), ("B", SEG_B), ("A", SEG_B)]);
    assert_error(load(&path), "Segment 'A' is defined more than once");
}

/// A segments file with no records at all cannot classify anything, so say so.
#[test]
fn empty_segments_file_is_rejected() {
    let dir = TempDir::new().unwrap();
    let path = dir.path().join("refs.fasta");
    File::create(&path).unwrap();
    assert_error(load(&path), "No segments found");
}

/// Sequences are compared as raw bytes, so a lower-case segment would once never match
/// an upper-case read and the tool silently reported nothing. Both sides are now
/// upper-cased on load, and a fully lower-case FASTA and FASTQ classify normally.
#[test]
fn lowercase_segments_and_reads_are_handled() {
    let dir = TempDir::new().unwrap();
    let refs = dir.path().join("refs.fasta");
    let reads = dir.path().join("reads.fastq");
    let read = format!("{}{}{}", SEG_A, SPACER, SEG_B);
    write_fasta(&refs, &[("A", &SEG_A.to_lowercase()), ("B", SEG_B)]);
    write_fastq(&reads, &[("r", &read.to_lowercase())]);

    let segments = load(&refs).unwrap();
    assert_eq!(
        segments["A"].seq,
        SEG_A.as_bytes(),
        "segments are upper-cased"
    );
    let loaded = load_reads(&reads).unwrap().reads;
    assert_eq!(as_pairs(&loaded), vec![("r".to_string(), read.clone())]);

    let (classified, _) = process_reads(loaded, &segments, (0, 0), false, false, false).unwrap();
    assert_eq!(classified[0].segments, "A-B");
}

/// An aligned BAM repeats a read once per alignment and hard-clips supplementary
/// records, so counts and sequences would both be wrong. Mapped records are refused
/// with an error that says what to do about it.
#[test]
fn aligned_bam_is_rejected() {
    let dir = TempDir::new().unwrap();
    let path = dir.path().join("aligned.bam");
    // No UNMAPPED flag means the record is mapped to a reference.
    write_bam(&path, &[("r1", SEG_A, Flags::empty())]);
    assert_error(load_reads(&path), "is an aligned BAM file");
    assert_error(load_reads(&path), "samtools fastq");
}

/// Secondary and supplementary alignments are the specific way an aligned BAM double
/// counts reads, so they are refused even alongside the unmapped flag.
#[test]
fn secondary_and_supplementary_bam_records_are_rejected() {
    let dir = TempDir::new().unwrap();
    for (name, flag) in [
        ("secondary.bam", Flags::SECONDARY),
        ("supplementary.bam", Flags::SUPPLEMENTARY),
    ] {
        let path = dir.path().join(name);
        write_bam(&path, &[("r1", SEG_A, Flags::UNMAPPED | flag)]);
        assert_error(load_reads(&path), "is an aligned BAM file");
    }
}

/// A BAM record that is individually unusable is skipped with a warning naming it,
/// while the rest of the file still loads.
#[test]
fn unusable_bam_records_are_skipped_not_fatal() {
    let dir = TempDir::new().unwrap();
    let path = dir.path().join("reads.bam");
    write_bam(
        &path,
        &[
            ("good1", SEG_A, Flags::UNMAPPED),
            ("no_sequence", "", Flags::UNMAPPED),
            ("good2", SEG_B, Flags::UNMAPPED),
        ],
    );
    assert_eq!(
        as_pairs(&load_reads(&path).unwrap().reads),
        vec![
            ("good1".to_string(), SEG_A.to_string()),
            ("good2".to_string(), SEG_B.to_string()),
        ]
    );
}

/// A malformed FASTQ record must not cost the rest of the file. This input holds a
/// stray line between records and a truncated record at the end; both are skipped with
/// a warning and every well-formed read either side of them is still returned.
#[test]
fn malformed_fastq_records_are_skipped_not_fatal() {
    let dir = TempDir::new().unwrap();
    let path = dir.path().join("reads.fastq");
    let mut f = File::create(&path).unwrap();
    write!(
        f,
        "@good1\n{SEG_A}\n+\nIIIIIIIIIIIIIIIIIIII\n\
         this stray line is not a record header\n\
         @good2\n{SEG_B}\n+\nIIIIIIIIIIIIIIIIIIII\n\
         @truncated\n{SEG_A}\n"
    )
    .unwrap();
    drop(f);

    assert_eq!(
        as_pairs(&load_reads(&path).unwrap().reads),
        vec![
            ("good1".to_string(), SEG_A.to_string()),
            ("good2".to_string(), SEG_B.to_string()),
        ]
    );
}

/// Asking for anchoring without supplying the anchors is a setup mistake, and the
/// message says which segment is missing and how to proceed.
#[test]
fn missing_anchor_segments_are_reported_clearly() {
    let segments = seg_map(&[("A", SEG_A)]);
    let reads = vec![raw("r", &construct())];
    let result = process_reads(reads, &segments, (0, 0), false, true, false);
    assert_error(result, "no segment named 'start'");
}

// ---------------------------------------------------------------------------
// Read excision and reorientation
// ---------------------------------------------------------------------------

const REORIENT_SEQ: &str = "ACGTACGGTT";

/// start < end: the region is simply excised and the orientation left alone.
#[test]
fn excises_forward_region() {
    assert_eq!(cut(REORIENT_SEQ, 2, 6, false), "GTAC");
}

/// start > end with reverse complementing requested: the read is reverse complemented
/// and the coordinates mirrored. This is how a reverse-strand read, whose start anchor
/// sits to the right of its end anchor, is restored to construct orientation.
#[test]
fn excises_and_reverse_complements_backwards_region() {
    assert_eq!(cut(REORIENT_SEQ, 8, 3, true), "CCGTA");
}

/// start > end without reverse complementing: the read is only reversed, not
/// complemented, so the two backwards paths stay distinguishable.
#[test]
fn excises_backwards_region_without_complementing() {
    assert_eq!(cut(REORIENT_SEQ, 8, 3, false), "GGCAT");
}

// ---------------------------------------------------------------------------
// Alignment scoring
// ---------------------------------------------------------------------------
//
// These pin the numbers every later decision is made from: a match is +2, a mismatch
// -1, a gap -2, and the normalised score is the raw score divided by segment length.
// The score threshold passed to each call is set to the expected peak so that only the
// intended alignment qualifies, which also pins the score exactly.

/// A perfect 20 bp match scores +2 per base and normalises to the 2.0 maximum, and the
/// reported span is exactly the region of the read the segment covers.
#[test]
fn scores_a_perfect_match() {
    let read = format!("{}{}{}", SPACER, SEG_A, SPACER);
    assert_eq!(only(align(SEG_A, &read, 40)), (40, 2.0, 12, 32));
}

/// One substitution costs 3 relative to a match (+2 becomes -1), so a 20 bp segment
/// with a single mismatch scores 37 across the same 20 bp span.
#[test]
fn scores_a_mismatch() {
    let read = format!("{}{}{}", SPACER, mutate(SEG_A, &[10]), SPACER);
    assert_eq!(only(align(SEG_A, &read, 37)), (37, 37.0 / 20.0, 12, 32));
}

/// A base missing from the read forces a gap: 19 matches (+38) less one gap (-2) = 36,
/// covering only 19 bases of the read.
#[test]
fn scores_a_deletion_in_the_read() {
    let read = format!("{}{}{}{}", SPACER, &SEG_A[..10], &SEG_A[11..], SPACER);
    assert_eq!(only(align(SEG_A, &read, 36)), (36, 36.0 / 20.0, 12, 31));
}

/// A base inserted into the read also costs one gap: 20 matches (+40) less 2 = 38,
/// stretched across 21 bases of the read.
#[test]
fn scores_an_insertion_in_the_read() {
    let read = format!("{}{}G{}{}", SPACER, &SEG_A[..10], &SEG_A[10..], SPACER);
    assert_eq!(only(align(SEG_A, &read, 38)), (38, 38.0 / 20.0, 12, 33));
}

/// Alignment is semi-global: gaps at either end of the read are free, so a segment
/// flush against the start or the end of a read scores the same as one in the middle.
#[test]
fn scores_segments_flush_with_either_end_of_the_read() {
    let at_start = format!("{}{}", SEG_A, SPACER);
    let at_end = format!("{}{}", SPACER, SEG_A);
    assert_eq!(only(align(SEG_A, &at_start, 40)), (40, 2.0, 0, 20));
    assert_eq!(only(align(SEG_A, &at_end, 40)), (40, 2.0, 12, 32));
}

/// The score threshold is inclusive: an alignment scoring exactly the minimum is kept,
/// and raising the minimum by one discards it.
#[test]
fn score_threshold_is_inclusive() {
    let read = format!("{}{}{}", SPACER, mutate(SEG_A, &[10]), SPACER);
    assert_eq!(align(SEG_A, &read, 37).len(), 1);
    assert!(align(SEG_A, &read, 38).is_empty());
}

/// Every occurrence of a repeated segment is reported separately, in read order. This
/// is the point of the multiple-traceback aligner: one segment, many hits per read.
#[test]
fn reports_every_occurrence_of_a_repeated_segment() {
    let read = format!("{}{}{}", SEG_A, SPACER, SEG_A);
    let spans: Vec<_> = align(SEG_A, &read, 40)
        .iter()
        .map(|a| (a.target_start, a.target_end))
        .collect();
    assert_eq!(spans, vec![(0, 20), (32, 52)]);
}

/// The separate aligner used to locate the start/end anchors scores identically and
/// reports the read span the anchor covers, which is what trimming is driven from.
#[test]
fn scores_and_locates_the_start_anchor() {
    let read = format!("{}{}{}", JUNK_5, START, SEG_A);
    let alignment = best_semiglobal(START.as_bytes(), read.as_bytes(), ANCHOR_SCORING);
    assert_eq!(alignment.score, 40);
    assert_eq!((alignment.ystart, alignment.yend), (10, 30));

    let degraded = format!("{}{}{}", JUNK_5, mutate(START, &[5]), SEG_A);
    let alignment = best_semiglobal(START.as_bytes(), degraded.as_bytes(), ANCHOR_SCORING);
    assert_eq!(alignment.score, 37);
}

// ---------------------------------------------------------------------------
// Segment classification and overlap resolution
// ---------------------------------------------------------------------------

/// With no anchors supplied and nothing found, the segment string is empty - the read
/// is still reported, but nothing is claimed about it.
#[test]
fn read_without_segments_reports_nothing() {
    let read = format!("{0}{0}{0}", SPACER);
    assert_eq!(classify(&[("A", SEG_A), ("B", SEG_B)], &read), "");
}

/// Well-separated segments are all reported, ordered by their position in the read.
/// B is planted first here, so a passing result also proves the ordering is positional
/// rather than following the order the segments were defined in.
#[test]
fn separated_segments_are_reported_in_read_order() {
    let read = format!("{}{}{}", SEG_B, SPACER, SEG_A);
    assert_eq!(classify(&[("A", SEG_A), ("B", SEG_B)], &read), "B-A");
}

/// A segment found on the opposite strand is reported under a starred name, so the
/// orientation of each individual segment survives into the output.
#[test]
fn reverse_complement_segments_are_starred() {
    let read = format!("{}{}{}", SEG_A, SPACER, rc(SEG_B));
    assert_eq!(classify(&[("A", SEG_A), ("B", SEG_B)], &read), "A-B*");
}

/// A segment whose best alignment falls below the score threshold is discarded. The
/// planted copy of B carries five mismatches, scoring 25 against a threshold of 30.
#[test]
fn segments_below_the_score_threshold_are_discarded() {
    let degraded = mutate(SEG_B, &[2, 6, 10, 14, 18]);
    let read = format!("{}{}{}", SEG_A, SPACER, degraded);
    assert_eq!(classify(&[("A", SEG_A), ("B", SEG_B)], &read), "A");
}

/// The 'start'/'end' bounds appear in the segment string only when the user supplied
/// those anchors. On the identical read and segment set, classifying as-given reports
/// just what was found, while classifying an anchor-trimmed read reports the bounds.
#[test]
fn anchors_bound_the_segment_string_only_when_supplied() {
    let read = format!("{}{}{}", SEG_A, SPACER, SEG_B);
    assert_eq!(classify(&[("A", SEG_A), ("B", SEG_B)], &read), "A-B");
    assert_eq!(
        classify_anchored(&[("A", SEG_A), ("B", SEG_B)], &read),
        "start-A-B-end"
    );
}

/// Without anchors, segments the user happens to have named 'start' and 'end' are not
/// special: they are searched for and reported like any other segment, rather than
/// silently ignored.
#[test]
fn segments_named_start_and_end_are_ordinary_when_not_anchoring() {
    let read = format!("{}{}{}", START, SPACER, SEG_A);
    assert_eq!(
        classify(&[("start", START), ("end", END), ("A", SEG_A)], &read),
        "start-A"
    );
}

/// With anchors supplied the read has already been trimmed to run from 'start' to
/// 'end', so those two are reported as the bounds and never again from within the
/// read - even though the trimmed read still physically contains both sequences.
#[test]
fn anchor_segments_are_not_reported_within_an_anchored_read() {
    let read = format!("{}{}{}", START, SEG_A, END);
    assert_eq!(
        classify_anchored(&[("start", START), ("end", END), ("A", SEG_A)], &read),
        "start-A-end"
    );
}

/// The per-segment threshold is min_norm_score x segment length, so at 1.5 a 20 bp
/// segment must score 30. Three mismatches (score 31) are tolerated; four (28) are not.
#[test]
fn classification_threshold_tracks_the_normalised_score() {
    let three = format!("{}{}", SPACER, mutate(SEG_B, &[2, 8, 14]));
    assert_eq!(classify(&[("B", SEG_B)], &three), "B");

    let four = format!("{}{}", SPACER, mutate(SEG_B, &[2, 8, 14, 18]));
    assert_eq!(classify(&[("B", SEG_B)], &four), "");
}

/// A segment occurring twice in one read is reported twice: the near-identical
/// alignments around each peak collapse, but the two genuine hits both survive.
#[test]
fn repeated_segments_are_reported_once_per_occurrence() {
    let read = format!("{}{}{}", SEG_A, SPACER, SEG_A);
    assert_eq!(classify(&[("A", SEG_A)], &read), "A-A");
}

/// Assert `segment` aligns across exactly `span` in `read` scoring exactly `score`, and
/// that this clears the classifier's threshold. Used by the repeated-segment tests below
/// to establish that every planted copy is a genuine above-threshold hit, so whatever
/// the classifier drops is dropped by the overlap rules and not by the score cutoff.
fn assert_hit_above_threshold(segment: &str, read: &str, span: (usize, usize), score: i32) {
    let threshold = (segment.len() as f32 * MIN_SCORE) as i32;
    assert!(
        score >= threshold,
        "the copy at {span:?} scores {score}, which is below the threshold of {threshold}"
    );
    let found = align(segment, read, threshold);
    assert!(
        found
            .iter()
            .any(|a| (a.target_start, a.target_end) == span && a.score == score),
        "expected a hit at {span:?} scoring {score}, but found {:?}",
        found
            .iter()
            .map(|a| (a.target_start, a.target_end, a.score))
            .collect::<Vec<_>>()
    );
}

/// Three well-separated copies of one segment are all reported. Each is a perfect match,
/// so all three clear the threshold and none overlaps another.
#[test]
fn repeated_segment_not_overlapping_keeps_every_copy() {
    let read = format!("{0}{1}{0}{1}{0}", SEG_A, SPACER);
    assert_hit_above_threshold(SEG_A, &read, (0, 20), 40);
    assert_hit_above_threshold(SEG_A, &read, (32, 52), 40);
    assert_hit_above_threshold(SEG_A, &read, (64, 84), 40);
    assert_eq!(classify(&[("A", SEG_A)], &read), "A-A-A");
}

/// Two copies of one segment overlapping by exactly 3 bp - the tolerance - are both
/// kept. Appending the segment from offset 3 makes the second copy start 3 bases before
/// the first one ends; it carries two mismatches where the copies interlock but still
/// scores 34 against a threshold of 30, so both are genuine hits.
#[test]
fn repeated_segment_overlapping_within_tolerance_keeps_both_copies() {
    let read = format!("{}{}", SEG_A, &SEG_A[3..]);
    assert_hit_above_threshold(SEG_A, &read, (0, 20), 40);
    assert_hit_above_threshold(SEG_A, &read, (17, 37), 34);
    assert_eq!(classify(&[("A", SEG_A)], &read), "A-A");
}

/// Copies overlapping by more than the tolerance cannot both be real, so the
/// better-scoring one survives. This needs SEG_REPEAT rather than SEG_A: where two
/// copies overlap they must agree on the bases they share, and for a non-repetitive
/// segment that is impossible - each extra base of overlap costs the second copy roughly
/// a mismatch, so it drops below the threshold before the overlap gets interesting.
/// SEG_REPEAT begins and ends with the same six bases, so copies placed six apart agree
/// across the join. Mismatches are then introduced only in the tail, which lies beyond
/// the first copy, so it scores 40 and the second 34.
#[test]
fn repeated_segment_overlapping_beyond_tolerance_keeps_the_best_copy() {
    let read = format!("{}{}", SEG_REPEAT, mutate(&SEG_REPEAT[6..], &[2, 9]));
    assert_hit_above_threshold(SEG_REPEAT, &read, (0, 20), 40);
    assert_hit_above_threshold(SEG_REPEAT, &read, (14, 34), 34);
    assert_eq!(classify(&[("A", SEG_REPEAT)], &read), "A");
}

/// One read exercising both outcomes at once: a perfect copy, a second overlapping it by
/// 6 bp that loses on score, and a third clear of both that is kept. Two of the three
/// above-threshold copies survive.
#[test]
fn repeated_segment_mixes_dropped_and_kept_copies_in_one_read() {
    let read = format!(
        "{}{}{}{}",
        SEG_REPEAT,
        mutate(&SEG_REPEAT[6..], &[2, 9]),
        SPACER,
        SEG_REPEAT
    );
    assert_hit_above_threshold(SEG_REPEAT, &read, (0, 20), 40);
    assert_hit_above_threshold(SEG_REPEAT, &read, (14, 34), 34);
    assert_hit_above_threshold(SEG_REPEAT, &read, (46, 66), 40);
    assert_eq!(classify(&[("A", SEG_REPEAT)], &read), "A-A");
}

/// Two segments that overlap by 3 bp - the maximum tolerated - are both kept. SEG_A
/// ends with "CAT" and SEG_D is built to begin with it, so their alignments share
/// exactly three bases of the read and neither should suppress the other.
#[test]
fn segments_overlapping_within_tolerance_are_both_kept() {
    let seg_d = format!("CAT{}", &SEG_B[..17]);
    let read = format!("{}{}", SEG_A, &seg_d[3..]);
    assert_eq!(
        classify(&[("A", SEG_A), ("D", seg_d.as_str())], &read),
        "A-D"
    );
}

/// When two segments overlap by more than the 3 bp tolerance only one can survive, and
/// it must be the better-scoring alignment. SEG_A matches perfectly here while SEG_E,
/// which covers a window overlapping it by 12 bp, carries two mismatches.
#[test]
fn heavily_overlapping_segments_resolve_to_the_better_alignment() {
    let read = format!("{}GGTTAACC", SEG_A);
    let seg_e = mutate(&read[8..28], &[3, 15]);
    assert_eq!(classify(&[("A", SEG_A), ("E", seg_e.as_str())], &read), "A");
}

/// The mirror of the previous case on the identical read: with the two mismatches moved
/// onto SEG_A instead, SEG_E now scores higher and is the one retained. Together the
/// pair show the survivor is picked by alignment quality, not by position or ordering.
#[test]
fn heavily_overlapping_segments_resolve_by_score_not_position() {
    let read = format!("{}GGTTAACC", SEG_A);
    let seg_e = read[8..28].to_string();
    let seg_a_degraded = mutate(SEG_A, &[4, 15]);
    assert_eq!(
        classify(
            &[("A", seg_a_degraded.as_str()), ("E", seg_e.as_str())],
            &read
        ),
        "E"
    );
}

/// A long segment that scores above threshold but completely contains a shorter,
/// perfectly matching one is discarded in favour of the contained alignment. This
/// exercises the containment rule specifically, rather than the partial-overlap rule.
#[test]
fn containing_alignment_is_discarded_in_favour_of_better_contained_one() {
    let short = "GCTACGGATCATGG";
    let read = format!("CGATGC{}TTAACC", short);
    let long = mutate(&read, &[2, 12, 22]);
    assert_eq!(
        classify(&[("SHORT", short), ("LONG", long.as_str())], &read),
        "SHORT"
    );
}

// ---------------------------------------------------------------------------
// Per-segment score thresholds
// ---------------------------------------------------------------------------

/// Write a segments file and load it with --per-segment-scores on.
fn load_scored_fasta(dir: &TempDir, records: &[(&str, &str)]) -> Result<HashMap<String, Segment>> {
    let path = dir.path().join("refs.fasta");
    write_fasta(&path, records);
    load_scored(&path)
}

/// The threshold in square brackets becomes that segment's own, and the brackets are not
/// part of the name it is reported under. A segment without them takes the global value,
/// so the two forms can be mixed in one file.
#[test]
fn a_bracketed_score_sets_that_segments_threshold_and_leaves_the_name_clean() {
    let dir = TempDir::new().unwrap();
    let segments = load_scored_fasta(&dir, &[("A[1.9]", SEG_A), ("B", SEG_B)]).unwrap();

    assert_eq!(segments.len(), 2);
    assert_eq!(segments["A"].seq, SEG_A.as_bytes());
    assert_eq!(segments["A"].min_norm_score, 1.9);
    // B said nothing, so it is held to whatever --min-norm-score was.
    assert_eq!(segments["B"].min_norm_score, MIN_SCORE);
}

/// A bracketed value always comes off the name, so a segment is reported under the same
/// name however the run is configured. Without the flag the value inside is simply not
/// used, and the segment falls back to the global threshold.
#[test]
fn brackets_come_off_the_name_even_without_the_flag() {
    let dir = TempDir::new().unwrap();
    let path = dir.path().join("refs.fasta");
    write_fasta(&path, &[("A[1.9]", SEG_A)]);

    let segments = load(&path).unwrap();
    assert_eq!(segments["A"].seq, SEG_A.as_bytes());
    assert_eq!(
        segments["A"].min_norm_score, MIN_SCORE,
        "1.9 is not applied"
    );
    assert!(!segments.contains_key("A[1.9]"));

    // The same file with the flag on: same name, and now the score is used.
    assert_eq!(load_scored(&path).unwrap()["A"].min_norm_score, 1.9);
}

/// Without the flag nothing inside the brackets is interpreted, so a bracketed name that
/// is not a score is dropped rather than rejected. It only has to be a valid score when
/// the flag says it is going to be read as one.
#[test]
fn a_bracketed_value_is_not_validated_unless_it_is_going_to_be_used() {
    let dir = TempDir::new().unwrap();
    let path = dir.path().join("refs.fasta");
    write_fasta(&path, &[("gene[human]", SEG_A)]);

    assert_eq!(load(&path).unwrap()["gene"].seq, SEG_A.as_bytes());
    assert_error(load_scored(&path), "'human' is not a number");
}

/// Because the brackets come off first, two records differing only inside them are one
/// segment named twice - with or without the flag.
#[test]
fn records_differing_only_inside_the_brackets_are_a_duplicate() {
    let dir = TempDir::new().unwrap();
    let path = dir.path().join("refs.fasta");
    write_fasta(&path, &[("A[1.5]", SEG_A), ("A[1.7]", SEG_B)]);
    assert_error(load(&path), "Segment 'A' is defined more than once");
    assert_error(load_scored(&path), "Segment 'A' is defined more than once");
}

/// The threshold reaches classification, which is the whole point. At 20 bp a single
/// mismatch scores 37, or 1.85 normalised: above 1.8 and below 1.9. One read and one
/// sequence therefore classify differently depending only on the bracketed score.
#[test]
fn a_bracketed_score_decides_whether_a_segment_is_found() {
    let read = format!("{}{}{}", SPACER, mutate(SEG_A, &[10]), SPACER);
    assert_hit_above_threshold(SEG_A, &read, (12, 32), 37);

    let classify_at = |min_norm_score: f32| {
        segment_string(&classify_read_segments(
            &seg_map_scored(&[("A", SEG_A, min_norm_score)]),
            read.as_bytes(),
            false,
        ))
    };
    assert_eq!(classify_at(1.8), "A", "37/20 = 1.85 clears 1.8");
    assert_eq!(classify_at(1.9), "", "...and falls short of 1.9");
}

/// Segments in one file are found at their own thresholds independently, so a strict
/// segment and a lenient one can sit side by side in the same read. Both carry one
/// mismatch and so both score 1.85; only the lenient one is reported.
#[test]
fn each_segment_is_found_at_its_own_threshold() {
    let read = format!(
        "{}{}{}{}",
        SPACER,
        mutate(SEG_A, &[10]),
        mutate(SEG_B, &[10]),
        SPACER
    );
    let segments = seg_map_scored(&[("strict", SEG_A, 1.9), ("lenient", SEG_B, 1.8)]);
    assert_eq!(
        segment_string(&classify_read_segments(&segments, read.as_bytes(), false)),
        "lenient"
    );
}

/// Anchors are segments too, so a bracketed score applies to them as well - and to each
/// anchor separately. A read whose 'start' carries one mismatch is dropped when 'start'
/// demands 1.9 and kept when it demands 1.8, while 'end' is untouched either way.
#[test]
fn a_bracketed_score_applies_to_the_anchors() {
    let degraded = format!(
        "{}{}{}{}{}{}",
        JUNK_5,
        mutate(START, &[10]),
        SEG_A,
        SPACER,
        END,
        JUNK_3
    );
    let anchored_at = |start_score: f32| {
        process_reads(
            vec![raw("read", &degraded)],
            &seg_map_scored(&[
                ("start", START, start_score),
                ("end", END, MIN_SCORE),
                ("A", SEG_A, MIN_SCORE),
            ]),
            (0, 0),
            false,
            true,
            false,
        )
        .unwrap()
    };
    let (kept, _) = anchored_at(1.8);
    assert_eq!(kept.len(), 1, "a degraded start clears 1.8");
    assert_eq!(kept[0].segments, "start-A-end");

    let (dropped, summary) = anchored_at(1.9);
    assert!(dropped.is_empty(), "...and falls short of 1.9");
    assert_eq!(
        summary.start_anchor_not_found, 1,
        "the start is what failed"
    );
}

/// Brackets that cannot be read as a score are a mistake worth stopping for. Quietly
/// treating them as part of the name would leave a segment that is never found and no
/// indication of why.
///
/// Each way of getting it wrong gets its own message, because "invalid score" tells you
/// nothing you did not already know. The table is the specification: the left column is
/// what a user might write and the right is the phrase that has to come back.
#[test]
fn every_way_of_writing_a_bad_score_gets_its_own_message() {
    let dir = TempDir::new().unwrap();
    let cases = [
        ("A[high]", "'high' is not a number"),
        ("A[1.5.2]", "'1.5.2' is not a number"),
        ("A[]", "the square brackets are empty"),
        // NaN parses as a float but is neither too high nor too low, so it is turned away
        // before the range check, whose advice would make no sense for it.
        ("A[NaN]", "'NaN' is not a usable score"),
        ("A[inf]", "'inf' is not a usable score"),
        // Out of range says what would have happened, not just what the bounds are.
        ("A[2.5]", "no alignment can score above 2.0"),
        ("A[-0.5]", "the segment would be found everywhere"),
        // A score with nothing in front of it names no segment.
        ("[1.5]", "there is no name in front of the brackets"),
        // Brackets that were opened and never closed are a typo, not a name.
        ("A[1.5", "the square bracket is never closed"),
        (
            "A[1.5]x",
            "the score has to come last, but 'x' follows the closing bracket",
        ),
    ];
    for (id, expected) in cases {
        assert_error(load_scored_fasta(&dir, &[(id, SEG_A)]), expected);
    }
    // A FASTA ID is the first word of the header, so an ID with spaces inside the
    // brackets cannot reach here through a file. The parser still has to handle it.
    assert_error(
        parse_segment_name("A[   ]", true),
        "the square brackets are empty",
    );
}

/// Whatever went wrong, the message says which record it was: the raw ID as written, the
/// file, and the record number. With a malformed bracket there may be no usable name to
/// report, and in a file of hundreds the number is what makes it findable.
#[test]
fn a_bad_score_is_reported_against_the_record_it_came_from() {
    let dir = TempDir::new().unwrap();
    let error = load_scored_fasta(&dir, &[("ok", SEG_A), ("B[oops]", SEG_B)]).unwrap_err();
    assert!(error.contains("Segment 'B[oops]'"), "{error}");
    assert!(error.contains("refs.fasta"), "{error}");
    assert!(error.contains("(record 2)"), "{error}");
}

/// The range is stated as the same 0.0 to 2.0 the global --min-norm-score is held to, so
/// the two cannot drift apart in a user's head.
#[test]
fn a_bad_score_says_how_to_write_a_good_one() {
    let dir = TempDir::new().unwrap();
    let error = load_scored_fasta(&dir, &[("A[high]", SEG_A)]).unwrap_err();
    assert!(error.contains("--per-segment-scores"), "{error}");
    assert!(error.contains("between 0.0 and 2.0"), "{error}");
    assert!(error.contains("'SEG1[1.5]'"), "names the syntax: {error}");
}

/// Scores are reported back as written rather than as the reparsed float, so the message
/// quotes what is actually in the file.
#[test]
fn a_bad_score_is_quoted_as_it_was_written() {
    let dir = TempDir::new().unwrap();
    let error = load_scored_fasta(&dir, &[("A[2.50]", SEG_A)]).unwrap_err();
    assert!(error.contains("'2.50'"), "not reformatted to 2.5: {error}");
}

/// Writing the syntax and forgetting the flag leaves a segment literally named
/// 'SEG1[1.5]' that never matches what was meant, with nothing in the run to say why. It
/// is not an error - the name is legal - but it is worth pointing at.
#[test]
fn a_score_written_without_the_flag_is_recognised_as_such() {
    assert!(looks_like_a_scored_name("SEG1[1.5]"));
    assert!(looks_like_a_scored_name("SEG1[0]"));
    // A bracketed name that is not a number is just a name, and must not be flagged.
    assert!(!looks_like_a_scored_name("gene[human]"));
    assert!(!looks_like_a_scored_name("SEG1[]"));
    assert!(!looks_like_a_scored_name("[1.5]"), "no name in front of it");
    assert!(!looks_like_a_scored_name("SEG1"));
    assert!(!looks_like_a_scored_name("SEG1[1.5"), "never closed");
}

/// The warning does not stop the file loading: the segment is still there, named without
/// its brackets as always, and held to the global threshold because the score it asked
/// for was not read.
#[test]
fn a_score_written_without_the_flag_still_loads() {
    let dir = TempDir::new().unwrap();
    let path = dir.path().join("refs.fasta");
    write_fasta(&path, &[("A[1.7]", SEG_A)]);
    let segments = load(&path).unwrap();
    assert_eq!(segments["A"].min_norm_score, MIN_SCORE);
}

/// A name ending in ']' with no '[' anywhere is not a bracket group, so it stays a name
/// whole. Only a matched pair at the end is taken off.
#[test]
fn only_a_closed_trailing_bracket_group_is_taken_off_the_name() {
    for flag in [false, true] {
        assert_eq!(
            parse_segment_name("odd]name]", flag).unwrap(),
            ("odd]name]".to_string(), None),
            "with the flag {flag}"
        );
    }
    // The last bracket pair wins, so a name may itself contain brackets.
    assert_eq!(
        parse_segment_name("SEG[1][1.5]", true).unwrap(),
        ("SEG[1]".to_string(), Some(1.5))
    );
    // Whitespace inside the brackets is tolerated rather than being a parse failure.
    assert_eq!(
        parse_segment_name("SEG[ 1.5 ]", true).unwrap(),
        ("SEG".to_string(), Some(1.5))
    );
    // Without the flag the same ID keeps its name and drops the value unread.
    assert_eq!(
        parse_segment_name("SEG[ 1.5 ]", false).unwrap(),
        ("SEG".to_string(), None)
    );
}

/// The ambiguity warning is judged against each segment's own threshold, so raising a
/// degenerate segment's score is a way to silence it - which is exactly what the warning
/// suggests doing.
#[test]
fn a_raised_threshold_settles_the_ambiguity_warning_for_that_segment() {
    let mostly_n = format!("{}{}", "N".repeat(16), &SEG_A[16..]);
    let flagged = |score: f32| -> Vec<String> {
        degenerate_segments(&seg_map_scored(&[
            ("vague", &mostly_n, score),
            ("specific", SEG_A, MIN_SCORE),
        ]))
        .into_iter()
        .map(|(name, _, _)| name)
        .collect()
    };
    // 16 Ns in 20 bases match 85% of random sequence: above the 83% that 1.5 demands...
    assert_eq!(flagged(MIN_SCORE), vec!["vague"]);
    // ...and below the 90% that 1.7 demands.
    assert!(flagged(1.7).is_empty());
}

// ---------------------------------------------------------------------------
// IUPAC ambiguity codes
// ---------------------------------------------------------------------------

// Bases match when the sets of nucleotides they name intersect, so a segment may use an
// ambiguity code to stand for several bases at once. On plain ACGT that rule is byte
// equality, which is what every test above this point relies on: those tests passing
// unchanged is what says nothing about unambiguous segments moved.

/// What each IUPAC code stands for, spelled out here independently of the table in
/// `main.rs` so the two have to agree rather than being the same statement twice.
const IUPAC_MEANINGS: [(u8, &str); 15] = [
    (b'A', "A"),
    (b'C', "C"),
    (b'G', "G"),
    (b'T', "T"),
    (b'R', "AG"),
    (b'Y', "CT"),
    (b'S', "CG"),
    (b'W', "AT"),
    (b'K', "GT"),
    (b'M', "AC"),
    (b'B', "CGT"),
    (b'D', "AGT"),
    (b'H', "ACT"),
    (b'V', "ACG"),
    (b'N', "ACGT"),
];

/// `SEG_A` with an `R` - A or G - where it has an A, so it matches `SEG_A` itself and
/// the variant carrying a G at that position, and nothing else.
const SEG_DEGENERATE: &str = "CGRTGCTAGCTACGGATCAT";
/// The G-carrying read sequence `SEG_DEGENERATE` matches but `SEG_A` does not.
const SEG_A_VARIANT: &str = "CGGTGCTAGCTACGGATCAT";

/// The match rule over every possible pair of bytes, checked against the meanings above.
/// Exhaustive because it is cheap: two bytes is 65,536 cases, and the interesting ones
/// are the edges - lower case, bytes that are not nucleotides at all, and each code
/// against each of the four bases.
#[test]
fn bases_match_exactly_when_the_codes_they_name_share_a_nucleotide() {
    let permits = |byte: u8| -> Option<&'static str> {
        IUPAC_MEANINGS
            .iter()
            .find(|(code, _)| *code == byte.to_ascii_uppercase())
            .map(|(_, bases)| *bases)
    };
    for query in 0..=u8::MAX {
        for target in 0..=u8::MAX {
            let want = match (permits(query), permits(target)) {
                (Some(a), Some(b)) => a.chars().any(|base| b.contains(base)),
                // Anything that is not a nucleotide code matches nothing at all, not
                // even a copy of itself.
                _ => false,
            };
            let got = BASE_MASK[query as usize] & BASE_MASK[target as usize] != 0;
            assert_eq!(
                got,
                want,
                "{} against {}",
                query.escape_ascii(),
                target.escape_ascii()
            );
        }
    }
}

/// An ambiguity code scores a full match against every base it permits and an ordinary
/// mismatch against the rest - exactly as if that base had been spelled out. Planting
/// each of the four bases in turn under an `R` shows both halves of that at once.
#[test]
fn an_ambiguity_code_scores_a_full_match_against_the_bases_it_permits() {
    let segment = format!("R{}", &SEG_A[1..]);
    for (base, score) in [('A', 40), ('G', 40), ('C', 37), ('T', 37)] {
        let read = format!("{}{}{}{}", SPACER, base, &SEG_A[1..], SPACER);
        assert_eq!(
            only(align(&segment, &read, score)),
            (score, score as f32 / 20.0, 12, 32),
            "read carrying {base} under the R"
        );
    }
}

/// The whole point of the feature: one degenerate segment finds both sequences it
/// stands for, and neither is favoured over the other.
#[test]
fn a_degenerate_segment_finds_every_sequence_it_stands_for() {
    for variant in [SEG_A, SEG_A_VARIANT] {
        let read = format!("{}{}{}", SPACER, variant, SPACER);
        assert_eq!(only(align(SEG_DEGENERATE, &read, 40)), (40, 2.0, 12, 32));
        assert_eq!(classify(&[("A", SEG_DEGENERATE)], &read), "A");
    }
}

/// A degenerate segment is still reported once per occurrence, and the two occurrences
/// need not be the same sequence. This is the repeat handling above meeting ambiguity
/// codes: nothing about overlap resolution changes, it just has more to resolve.
#[test]
fn a_degenerate_segment_is_reported_once_per_occurrence() {
    let read = format!("{}{}{}", SEG_A, SPACER, SEG_A_VARIANT);
    assert_hit_above_threshold(SEG_DEGENERATE, &read, (0, 20), 40);
    assert_hit_above_threshold(SEG_DEGENERATE, &read, (32, 52), 40);
    assert_eq!(classify(&[("A", SEG_DEGENERATE)], &read), "A-A");
}

/// Segments are searched on both strands, so a degenerate one has to complement
/// correctly: R (A or G) must become Y (C or T), not stay R and not fall back to N.
/// Both variants are found reversed, and starred to say which strand they sit on.
#[test]
fn a_degenerate_segment_is_found_on_the_reverse_strand() {
    for variant in [SEG_A, SEG_A_VARIANT] {
        let read = format!("{}{}{}", SPACER, rc(variant), SPACER);
        assert_eq!(classify(&[("A", SEG_DEGENERATE)], &read), "A*");
    }
    assert_eq!(rc(SEG_DEGENERATE), "ATGATCCGTAGCTAGCAYCG");
}

/// An ambiguity code stays a statement about which bases are allowed rather than a
/// licence to match anything: a base it excludes costs an ordinary mismatch, so the code
/// buys exactly one base's worth of freedom and no more. C and T are what R rules out.
#[test]
fn a_degenerate_segment_still_rejects_the_bases_it_excludes() {
    for excluded in ['C', 'T'] {
        let read = format!("{}CG{excluded}TGCTAGCTACGGATCAT{}", SPACER, SPACER);
        assert_eq!(
            only(align(SEG_DEGENERATE, &read, 37)),
            (37, 37.0 / 20.0, 12, 32),
            "read carrying {excluded} under the R"
        );
    }
    // That one mismatch is enough to decide a read that is otherwise marginal. Three
    // further mismatches leave a permitted base just above the threshold on 31, while an
    // excluded one falls just below it on 28 and is not reported at all.
    let read = |seq: &str| format!("{}{}{}", SPACER, mutate(seq, &[8, 12, 16]), SPACER);
    assert_eq!(
        classify(&[("A", SEG_DEGENERATE)], &read(SEG_A_VARIANT)),
        "A"
    );
    assert_eq!(
        classify(&[("A", SEG_DEGENERATE)], &read("CGCTGCTAGCTACGGATCAT")),
        ""
    );
}

/// Anchors go through a different aligner from the segments, so the rule has to hold
/// there too: a degenerate `start` trims and orients a read carrying either sequence it
/// permits, and the segments between the anchors are then reported as usual.
#[test]
fn a_degenerate_anchor_trims_and_orients_a_read() {
    // START with a W - A or T - where it has an A.
    let degenerate_start = "AGGCWTTCGAGCTTAACGGT";
    let segments = vec![
        ("start", degenerate_start),
        ("end", END),
        ("A", SEG_A),
        ("B", SEG_B),
    ];
    for variant in [START, "AGGCTTTCGAGCTTAACGGT"] {
        let body = format!("{}{}{}{}{}", variant, SEG_A, SPACER, SEG_B, END);
        let read = format!("{}{}{}", JUNK_5, body, JUNK_3);
        assert_eq!(
            run(vec![raw("read", &read)], &segments, false, true),
            expect("read", "start-A-B-end"),
            "anchor variant {variant}"
        );
    }
}

/// The rule is symmetric, so an N a basecaller left in a read matches whatever the
/// segment asks for there rather than costing a mismatch. This is a change: before
/// ambiguity codes were understood, an N in a read mismatched every segment base.
#[test]
fn an_ambiguity_code_in_a_read_matches_the_segment_base_it_permits() {
    let read = format!("{}N{}{}", SPACER, &SEG_A[1..], SPACER);
    assert_eq!(only(align(SEG_A, &read, 40)), (40, 2.0, 12, 32));
    assert_eq!(classify(&[("A", SEG_A)], &read), "A");
}

/// Reads are not validated the way segments are - noisy data is expected, and refusing a
/// whole file over one stray byte would be hostile - so a character that is not a
/// nucleotide has to cost that one position and nothing more. In particular it must not
/// match a copy of itself, which plain byte equality would have let it do.
#[test]
fn a_character_that_is_not_a_nucleotide_matches_nothing_in_a_read() {
    let read = format!("{}X{}{}", SPACER, &SEG_A[1..], SPACER);
    assert_eq!(only(align(SEG_A, &read, 37)), (37, 37.0 / 20.0, 12, 32));
    let segment = format!("X{}", &SEG_A[1..]);
    assert_eq!(only(align(&segment, &read, 37)), (37, 37.0 / 20.0, 12, 32));
}

/// Ambiguity is free under this scoring - a code matches at full score however many
/// bases it permits - so a segment can be degenerate enough to clear the threshold on
/// sequence chosen at random. That is a legitimate thing to search for, so it warns
/// rather than failing, but it must not pass unremarked.
///
/// At the default threshold a segment needs (1.5 + 1) / 3 = 83% of its bases to match,
/// while an N matches 100% of random bases and a spelled-out base 25%. Sixteen Ns in
/// twenty bases comes to 85% and is called out; fifteen comes to 81% and is not.
#[test]
fn segments_too_ambiguous_for_the_threshold_are_reported() {
    let mostly_n = format!("{}{}", "N".repeat(16), &SEG_A[16..]);
    let some_n = format!("{}{}", "N".repeat(15), &SEG_A[15..]);
    let flagged = |min_norm_score: f32| -> Vec<String> {
        let segments = seg_map_scored(&[
            ("A", SEG_A, min_norm_score),
            ("mostly_n", &mostly_n, min_norm_score),
            ("some_n", &some_n, min_norm_score),
        ]);
        degenerate_segments(&segments)
            .into_iter()
            .map(|(name, _, _)| name)
            .collect()
    };
    assert_eq!(flagged(MIN_SCORE), vec!["mostly_n"]);
    // The check is relative to the threshold, not absolute: drop the threshold far
    // enough and the milder segment matches random sequence often enough to qualify too.
    assert_eq!(flagged(0.5), vec!["mostly_n", "some_n"]);
    // ...and a strict enough threshold clears both.
    assert!(flagged(2.0).is_empty());
}

/// The warning reports the chance of matching random sequence, which is what makes it
/// actionable, so the number has to be right. A segment of nothing but N matches every
/// base, and one with no ambiguity at all matches a quarter of them.
#[test]
fn the_reported_chance_of_a_random_match_is_the_average_over_the_bases() {
    let all_n = "N".repeat(20);
    let half_n = format!("{}{}", "N".repeat(10), &SEG_A[10..]);
    assert_eq!(chance_match_fraction(all_n.as_bytes()), 1.0);
    assert_eq!(chance_match_fraction(SEG_A.as_bytes()), 0.25);
    // Ten bases matching everything and ten matching one base in four.
    assert_eq!(
        chance_match_fraction(half_n.as_bytes()),
        (10.0 + 10.0 * 0.25) / 20.0
    );
    // Each code in turn, against the bases it permits out of the four.
    for (code, bases) in IUPAC_MEANINGS {
        assert_eq!(
            chance_match_fraction(&[code]),
            bases.len() as f32 / 4.0,
            "code {}",
            code as char
        );
    }
}

// ---------------------------------------------------------------------------
// Anchoring and realignment on the start/end segments
// ---------------------------------------------------------------------------

/// A forward read carrying adapter junk on both flanks is trimmed back to START..END
/// before classification, so nothing outside the anchors can contribute segments.
#[test]
fn forward_read_is_trimmed_to_the_start_and_end_anchors() {
    let read = format!("{}{}{}", JUNK_5, construct(), JUNK_3);
    assert_eq!(
        run(vec![raw("r", &read)], &anchored_segments(), false, true),
        expect("r", "start-A-B-end")
    );
}

/// The reverse complement of that same read must be flipped back into construct
/// orientation, giving exactly the same segment string - not the reversed, starred
/// reading "start-B*-A*-end" that classifying it as-is would produce.
#[test]
fn reverse_read_is_reoriented_to_match_the_forward_result() {
    let read = rc(&format!("{}{}{}", JUNK_5, construct(), JUNK_3));
    assert_eq!(
        run(vec![raw("r", &read)], &anchored_segments(), false, true),
        expect("r", "start-A-B-end")
    );
}

/// A read lacking the START anchor cannot be oriented or trimmed reliably, so it is
/// dropped rather than reported from an arbitrary offset.
#[test]
fn read_missing_the_start_anchor_is_dropped() {
    let read = format!("{}{}{}{}{}{}", JUNK_5, SEG_A, SPACER, SEG_B, END, JUNK_3);
    assert!(run(vec![raw("r", &read)], &anchored_segments(), false, true).is_empty());
}

/// An anchor that is present but too degraded to clear the threshold also drops the
/// read, and does so on either strand - the forward and reverse branches each apply
/// the check. Six mismatches put END at 22 against a threshold of 30.
#[test]
fn read_with_a_degraded_anchor_is_dropped_on_either_strand() {
    let broken = format!(
        "{}{}{}{}{}",
        START,
        SEG_A,
        SPACER,
        SEG_B,
        mutate(END, &[1, 4, 7, 10, 13, 16])
    );
    let forward = format!("{}{}{}", JUNK_5, broken, JUNK_3);
    assert!(run(vec![raw("f", &forward)], &anchored_segments(), false, true).is_empty());
    assert!(
        run(
            vec![raw("r", &rc(&forward))],
            &anchored_segments(),
            false,
            true
        )
        .is_empty()
    );
}

// ---------------------------------------------------------------------------
// Circular constructs
// ---------------------------------------------------------------------------

/// A read off a circular template can begin anywhere, so END may precede START. This
/// rotation crosses the origin mid-construct: END-START-A-SPACER-B.
fn rotated_construct() -> String {
    format!("{}{}{}{}{}", END, START, SEG_A, SPACER, SEG_B)
}

/// Without --circular, a read whose START follows its END is inconsistent for a linear
/// construct and must be rejected rather than trimmed to a nonsensical span.
#[test]
fn wrapped_read_is_rejected_when_not_circular() {
    assert!(
        run(
            vec![raw("r", &rotated_construct())],
            &anchored_segments(),
            false,
            true
        )
        .is_empty()
    );
}

/// The reverse-strand mirror of the case above: a wrapped read sequenced off the
/// opposite strand is rejected without --circular too, so neither strand can slip
/// past the linear consistency check.
#[test]
fn wrapped_reverse_read_is_rejected_when_not_circular() {
    let read = rc(&rotated_construct());
    assert!(run(vec![raw("r", &read)], &anchored_segments(), false, true).is_empty());
}

/// Most reads off a circular template do not cross the origin. --circular must leave
/// those on the ordinary path rather than rotating them anyway, on either strand.
#[test]
fn unwrapped_reads_are_unaffected_by_the_circular_flag() {
    let forward = format!("{}{}{}", JUNK_5, construct(), JUNK_3);
    assert_eq!(
        run(vec![raw("r", &forward)], &anchored_segments(), true, true),
        expect("r", "start-A-B-end")
    );
    assert_eq!(
        run(
            vec![raw("r", &rc(&forward))],
            &anchored_segments(),
            true,
            true
        ),
        expect("r", "start-A-B-end")
    );
}

/// With --circular the very same read is un-rotated across the origin back into
/// START-A-B-END and classified normally.
#[test]
fn wrapped_forward_read_is_unrotated_when_circular() {
    assert_eq!(
        run(
            vec![raw("r", &rotated_construct())],
            &anchored_segments(),
            true,
            true
        ),
        expect("r", "start-A-B-end")
    );
}

/// A wrapped read sequenced off the opposite strand needs both un-rotation and reverse
/// complementing. It must land on the same answer as the forward wrapped read.
#[test]
fn wrapped_reverse_read_is_unrotated_and_reoriented_when_circular() {
    let read = rc(&rotated_construct());
    assert_eq!(
        run(vec![raw("r", &read)], &anchored_segments(), true, true),
        expect("r", "start-A-B-end")
    );
}

// ---------------------------------------------------------------------------
// Classification without start/end anchors
// ---------------------------------------------------------------------------

/// Without anchors the read is classified exactly as supplied, with no trimming, so
/// segments are still found even though flanking junk is left in place.
#[test]
fn unanchored_reads_are_classified_as_given() {
    let read = format!("{}{}{}{}", JUNK_5, SEG_A, SPACER, SEG_B);
    assert_eq!(
        run(
            vec![raw("r", &read)],
            &[("A", SEG_A), ("B", SEG_B)],
            false,
            false
        ),
        expect("r", "A-B")
    );
}

/// Without anchors there is nothing to orient against, so a reverse-strand read is
/// reported on the strand it arrived on: reversed order and starred names. This is the
/// contrast case for `reverse_read_is_reoriented_to_match_the_forward_result`.
#[test]
fn unanchored_reverse_reads_are_not_reoriented() {
    let read = rc(&format!("{}{}{}{}", JUNK_5, SEG_A, SPACER, SEG_B));
    assert_eq!(
        run(
            vec![raw("r", &read)],
            &[("A", SEG_A), ("B", SEG_B)],
            false,
            false
        ),
        expect("r", "B*-A*")
    );
}

// ---------------------------------------------------------------------------
// Read length filtering
// ---------------------------------------------------------------------------

/// Reads outside [min, max] are dropped before any alignment work, and the sentinel
/// bounds (0, 0) disable the check entirely rather than rejecting everything.
#[test]
fn reads_are_filtered_by_length() {
    let segments = seg_map(&[("A", SEG_A)]);
    let reads = || {
        vec![
            raw("short", SEG_A),                                  // 20 bp
            raw("ok", &format!("{}{}{}", SPACER, SEG_A, SPACER)), // 44 bp
            raw("long", &SEG_A.repeat(5)),                        // 100 bp
        ]
    };

    let bounded: Vec<_> = process_reads(reads(), &segments, (30, 50), false, false, false)
        .unwrap()
        .0
        .into_iter()
        .map(|r| r.name)
        .collect();
    assert_eq!(bounded, vec!["ok"]);

    let unbounded: Vec<_> = process_reads(reads(), &segments, (0, 0), false, false, false)
        .unwrap()
        .0
        .into_iter()
        .map(|r| r.name)
        .collect();
    assert_eq!(unbounded, vec!["short", "ok", "long"]);
}

// ---------------------------------------------------------------------------
// Equivalence of the linear-space aligner with the original
// ---------------------------------------------------------------------------
//
// `align_multiple` used to build the whole (query x target) matrix of scored cells,
// each recording which neighbours produced its score, then walk backwards from each
// end position to recover where the alignment began. It now carries the start position
// forward through two rolling rows instead. The tests below check that this is exactly
// equivalent and not merely close, by keeping the original implementation as a
// reference oracle and comparing every field of every alignment it produces.
//
// The oracle compares bases with `==` rather than by nucleotide set, so it only speaks
// for unambiguous input. That is deliberate - it is a copy of the original, not a copy
// of the current code - and it is why every generator below draws from A/T or ACGT. The
// IUPAC match rule is checked separately, against its own reference, in the section on
// ambiguity codes.

/// The original full-matrix implementation, preserved verbatim. Not used in production;
/// its only purpose is to be compared against.
mod reference {
    // Indexing style is preserved from the original implementation this mirrors;
    // rewriting it to satisfy the linter would defeat the point of a reference copy.
    #![allow(clippy::needless_range_loop)]

    use crate::MultiAlignment;

    #[derive(Clone, Copy)]
    struct Cell {
        score: i32,
        from_diagonal: bool,
        from_up: bool,
        from_left: bool,
    }

    impl Cell {
        /// A cell with no recorded predecessors yet.
        fn new(score: i32) -> Self {
            Cell {
                score,
                from_diagonal: false,
                from_up: false,
                from_left: false,
            }
        }
    }

    /// Fill the full scoring matrix, recording which neighbours produced each cell.
    fn build_dp_matrix(
        (match_score, mismatch_score, gap_score): (i32, i32, i32),
        query: &[u8],
        target: &[u8],
    ) -> Vec<Vec<Cell>> {
        let query_len = query.len();
        let target_len = target.len();
        let mut matrix = vec![vec![Cell::new(0); target_len + 1]; query_len + 1];
        for j in 1..=target_len {
            matrix[0][j] = Cell::new(0);
        }
        for i in 1..=query_len {
            matrix[i][0] = Cell::new(gap_score * i as i32);
            if i > 1 {
                matrix[i][0].from_up = true;
            }
        }
        for i in 1..=query_len {
            for j in 1..=target_len {
                let m = if query[i - 1] == target[j - 1] {
                    match_score
                } else {
                    mismatch_score
                };
                let diagonal_score = matrix[i - 1][j - 1].score + m;
                let up_score = matrix[i - 1][j].score + gap_score;
                let left_score = matrix[i][j - 1].score + gap_score;
                let max_score = diagonal_score.max(up_score).max(left_score);
                let mut cell = Cell::new(max_score);
                if diagonal_score == max_score {
                    cell.from_diagonal = true;
                }
                if up_score == max_score {
                    cell.from_up = true;
                }
                if left_score == max_score {
                    cell.from_left = true;
                }
                matrix[i][j] = cell;
            }
        }
        matrix
    }

    /// Score every end position, then trace each qualifying one back to its start.
    pub fn align_multiple(
        scoring: (i32, i32, i32),
        seg_name: &str,
        query: &[u8],
        target: &[u8],
        min_score: i32,
    ) -> Vec<MultiAlignment> {
        let query_len = query.len();
        let target_len = target.len();
        let matrix = build_dp_matrix(scoring, query, target);
        let mut alignments: Vec<MultiAlignment> = Vec::new();
        for end_j in 0..=target_len {
            let score = matrix[query_len][end_j].score;
            if score < min_score {
                continue;
            }
            // Greedy backward traceback preferring diagonal, then up, then left.
            let mut i = query_len;
            let mut j = end_j;
            while i > 0 && j > 0 {
                let cell = &matrix[i][j];
                if cell.from_diagonal {
                    i -= 1;
                    j -= 1;
                } else if cell.from_up {
                    i -= 1;
                } else if cell.from_left {
                    j -= 1;
                } else {
                    break;
                }
            }
            alignments.push(MultiAlignment {
                name: seg_name.to_string(),
                score,
                norm_score: score as f32 / query_len as f32,
                target_start: j,
                target_end: end_j,
                filter: false,
            });
        }
        alignments
    }
}

/// Scoring used throughout the comparison: the same +2 / -1 / -2 the classifier applies.
const SCORING: (i32, i32, i32) = (2, -1, -2);

/// A tiny deterministic PRNG, so the randomised comparison is reproducible and needs no
/// extra dependency.
struct Rng(u64);

impl Rng {
    /// xorshift64: cheap, deterministic, adequate for fixtures.
    fn next_u64(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.0 = x;
        x
    }

    /// A value in `0..n`.
    fn below(&mut self, n: usize) -> usize {
        (self.next_u64() % n as u64) as usize
    }

    /// A random sequence of `len` bases drawn from `alphabet`.
    fn seq(&mut self, len: usize, alphabet: &[u8]) -> Vec<u8> {
        (0..len)
            .map(|_| alphabet[self.below(alphabet.len())])
            .collect()
    }
}

/// Assert two alignment lists agree field for field.
fn assert_same_alignments(got: &[MultiAlignment], want: &[MultiAlignment], context: &str) {
    assert_eq!(
        got.len(),
        want.len(),
        "{context}: different number of alignments"
    );
    for (g, w) in got.iter().zip(want) {
        assert_eq!(g.name, w.name, "{context}: name");
        assert_eq!(g.score, w.score, "{context}: score");
        // Compared bitwise so the NaN from a zero-length query matches itself.
        assert_eq!(
            g.norm_score.to_bits(),
            w.norm_score.to_bits(),
            "{context}: norm_score"
        );
        assert_eq!(g.target_start, w.target_start, "{context}: target_start");
        assert_eq!(g.target_end, w.target_end, "{context}: target_end");
        assert_eq!(g.filter, w.filter, "{context}: filter");
    }
}

/// Exhaustive equivalence over every sequence pair up to 3 x 6 bases from a two-letter
/// alphabet. Short two-letter sequences are dense in ties - many cells have several
/// equally good predecessors - which is precisely where a tie-breaking difference
/// between the two implementations would surface as a different reported start.
#[test]
fn linear_space_aligner_matches_reference_exhaustively() {
    let aligner = MultiTracebackAligner::new(SCORING.0, SCORING.1, SCORING.2);
    let words = |len: usize| -> Vec<Vec<u8>> {
        (0..(1usize << len))
            .map(|mask| {
                (0..len)
                    .map(|bit| if mask >> bit & 1 == 1 { b'A' } else { b'T' })
                    .collect()
            })
            .collect()
    };
    let mut cases = 0usize;
    for query_len in 0..=3 {
        for query in words(query_len) {
            for target_len in 0..=6 {
                for target in words(target_len) {
                    for min_score in [-20, -4, -1, 0, 1, 2, 4, 6] {
                        let got = aligner.align_multiple("s", &query, &target, min_score);
                        let want =
                            reference::align_multiple(SCORING, "s", &query, &target, min_score);
                        assert_same_alignments(
                            &got,
                            &want,
                            &format!(
                                "query={:?} target={:?} min_score={min_score}",
                                String::from_utf8_lossy(&query),
                                String::from_utf8_lossy(&target)
                            ),
                        );
                        cases += 1;
                    }
                }
            }
        }
    }
    assert!(cases > 15_000, "expected a broad sweep, only ran {cases}");
}

/// Randomised equivalence at larger sizes over three alphabets: a single letter (every
/// position matches, so ties are everywhere), two letters, and realistic DNA.
#[test]
fn linear_space_aligner_matches_reference_on_random_inputs() {
    let aligner = MultiTracebackAligner::new(SCORING.0, SCORING.1, SCORING.2);
    let mut rng = Rng(0x9E3779B97F4A7C15);
    for alphabet in [&b"A"[..], &b"AT"[..], &b"ACGT"[..]] {
        for case in 0..400 {
            let query_len = rng.below(24);
            let target_len = rng.below(60);
            let query = rng.seq(query_len, alphabet);
            let target = rng.seq(target_len, alphabet);
            let min_score = rng.below(40) as i32 - 20;
            let got = aligner.align_multiple("s", &query, &target, min_score);
            let want = reference::align_multiple(SCORING, "s", &query, &target, min_score);
            assert_same_alignments(
                &got,
                &want,
                &format!(
                    "alphabet={} case={case} min_score={min_score}",
                    String::from_utf8_lossy(alphabet)
                ),
            );
        }
    }
}

/// Equivalence at a realistic scale: a 60 bp segment planted three times in a 2.4 kb
/// read - once cleanly, once with mismatches, once with a deleted base - so the
/// comparison covers gapped alignments and multiple hits, not just clean matches.
#[test]
fn linear_space_aligner_matches_reference_at_realistic_scale() {
    let aligner = MultiTracebackAligner::new(SCORING.0, SCORING.1, SCORING.2);
    let mut rng = Rng(0x0DDBA11);
    let segment = rng.seq(60, b"ACGT");

    let mut mismatched = segment.clone();
    mismatched[10] = if mismatched[10] == b'A' { b'C' } else { b'A' };
    mismatched[40] = if mismatched[40] == b'G' { b'T' } else { b'G' };
    let mut deleted = segment.clone();
    deleted.remove(25);

    let mut read = rng.seq(600, b"ACGT");
    for planted in [&segment, &mismatched, &deleted] {
        read.extend_from_slice(planted);
        read.extend(rng.seq(600, b"ACGT"));
    }

    for min_score in [-50, 0, 60, 90, 105, 114, 120] {
        let got = aligner.align_multiple("seg", &segment, &read, min_score);
        let want = reference::align_multiple(SCORING, "seg", &segment, &read, min_score);
        assert_same_alignments(&got, &want, &format!("min_score={min_score}"));
    }
}

/// Two equally good, mutually exclusive hits must resolve the same way on every run.
/// Segments live in a HashMap whose iteration order varies between runs, and the score
/// sort is stable, so without an explicit tiebreak the survivor would be whichever
/// happened to be visited first and the same input could give different answers. Both
/// segments here carry the same sequence, so they tie exactly and fully overlap; the
/// tiebreak keeps the alphabetically earlier name.
#[test]
fn tied_alignments_resolve_deterministically() {
    let read = format!("{}{}{}", SPACER, SEG_A, SPACER);
    let segments = [("A_first", SEG_A), ("Z_second", SEG_A)];
    assert_eq!(classify(&segments, &read), "A_first");
    // Each call builds a fresh HashMap, so a run-order dependence would show up here.
    for _ in 0..200 {
        assert_eq!(classify(&segments, &read), "A_first");
    }
}

// ---------------------------------------------------------------------------
// End-of-run summary
// ---------------------------------------------------------------------------

/// Every read is accounted for exactly once, under the reason that actually applies.
/// The eight reads here are built to land in the eight distinct outcomes, and the
/// classified and rejected counts must together add back up to the number of reads.
#[test]
fn summary_counts_every_read_under_the_right_reason() {
    let segments = seg_map(&anchored_segments());
    let reads = vec![
        raw("classified", &construct()),
        raw("tiny", "ACGT"),
        raw("huge", &SEG_A.repeat(20)),
        raw("no_start", &format!("{}{}", SEG_A, END)),
        raw("no_end", &format!("{}{}", START, SEG_A)),
        raw("neither", &SPACER.repeat(3)),
        raw(
            "opposite_strands",
            &format!("{}{}{}", START, SEG_A, rc(END)),
        ),
        raw("out_of_order", &rotated_construct()),
    ];
    let (classified, summary) =
        process_reads(reads, &segments, (10, 200), false, true, false).unwrap();

    assert_eq!(classified.len(), 1);
    assert_eq!(summary.classified, 1);
    assert_eq!(summary.too_short, 1, "4 bp read against a 10 bp minimum");
    assert_eq!(summary.too_long, 1, "400 bp read against a 200 bp maximum");
    assert_eq!(summary.start_anchor_not_found, 1);
    assert_eq!(summary.end_anchor_not_found, 1);
    assert_eq!(summary.neither_anchor_found, 1);
    assert_eq!(summary.anchors_on_opposite_strands, 1);
    assert_eq!(summary.anchors_out_of_order, 1);

    assert_eq!(summary.reads(), 8, "every read is accounted for");
    assert_eq!(summary.not_classified(), 7);
    assert_eq!(
        summary.classified + summary.not_classified(),
        summary.reads()
    );
}

/// Without anchoring, no anchor-related rejection can arise, so reads either classify or
/// fall outside the length bounds.
#[test]
fn summary_without_anchoring_only_reports_length_rejections() {
    let segments = seg_map(&[("A", SEG_A)]);
    let reads = vec![
        raw("ok", &format!("{}{}{}", SPACER, SEG_A, SPACER)),
        raw("tiny", "ACGT"),
    ];
    let (_, summary) = process_reads(reads, &segments, (10, 200), false, false, false).unwrap();
    assert_eq!(summary.classified, 1);
    assert_eq!(summary.too_short, 1);
    assert_eq!(summary.start_anchor_not_found, 0);
    assert_eq!(summary.end_anchor_not_found, 0);
    assert_eq!(summary.neither_anchor_found, 0);
    assert_eq!(summary.anchors_on_opposite_strands, 0);
    assert_eq!(summary.anchors_out_of_order, 0);
}

/// The rendered report lists only the categories that actually occurred, so a clean run
/// is not padded out with a column of zeros.
#[test]
fn summary_report_omits_categories_that_did_not_occur() {
    let summary = RunSummary {
        classified: 8,
        too_short: 2,
        ..Default::default()
    };
    let report = summary.render();
    assert!(report.contains("classified:"));
    assert!(report.contains("read too short"));
    assert!(
        !report.contains("read too long"),
        "zero counts are omitted:\n{report}"
    );
    assert!(
        !report.contains("out of order"),
        "zero counts are omitted:\n{report}"
    );
    assert!(report.contains("80.0%"), "percentages are shown:\n{report}");
}

// ---------------------------------------------------------------------------
// Classification counts CSV
// ---------------------------------------------------------------------------

/// Tally the given classifications and render the CSV --counts would write.
fn counts_csv(classifications: &[&str]) -> String {
    let mut counts = ClassificationCounts::default();
    for segments in classifications {
        counts.count(segments);
    }
    counts.render_csv()
}

/// Every classified read is counted exactly once under the string it produced, and the
/// rows come out most frequent first so the common result is the one you read.
#[test]
fn counts_csv_counts_each_classification_and_leads_with_the_commonest() {
    let csv = counts_csv(&["start-A-end", "start-B-end", "start-A-end", "start-A-end"]);
    assert_eq!(csv, "segments,count\nstart-A-end,3\nstart-B-end,1\n");
    // The counts account for every read handed in, so the file reconciles with the
    // number of lines in the results.
    let total: usize = csv
        .lines()
        .skip(1)
        .map(|line| line.rsplit_once(',').unwrap().1.parse::<usize>().unwrap())
        .sum();
    assert_eq!(total, 4);
}

/// Equally common classifications are ordered alphabetically rather than left to the
/// HashMap, whose iteration order varies between runs. Two runs over the same reads must
/// produce byte-identical files, or a diff between two conditions is unreadable.
#[test]
fn counts_csv_breaks_ties_alphabetically() {
    let classifications = ["start-C-end", "start-A-end", "start-B-end"];
    assert_eq!(
        counts_csv(&classifications),
        "segments,count\nstart-A-end,1\nstart-B-end,1\nstart-C-end,1\n"
    );
    // ...and shuffling the reads does not change the file.
    let mut reversed = classifications;
    reversed.reverse();
    assert_eq!(counts_csv(&reversed), counts_csv(&classifications));
}

/// A read that was classified but carried no recognisable segments is a real result -
/// quite different from a read that was dropped - so it gets a row of its own, with an
/// empty first field.
#[test]
fn counts_csv_keeps_reads_that_matched_no_segments() {
    assert_eq!(
        counts_csv(&["", "start-A-end", ""]),
        "segments,count\n,2\nstart-A-end,1\n"
    );
}

/// A run that classified nothing still writes a valid CSV rather than an empty file, so
/// whatever reads it stays a table with no rows instead of failing to parse.
#[test]
fn counts_csv_of_nothing_is_the_header_alone() {
    assert_eq!(counts_csv(&[]), "segments,count\n");
}

/// Segment names come from FASTA headers, so a comma or a quote in one is unlikely but
/// perfectly legal. Left alone it would shift every column after it, so fields are
/// quoted when they have to be and left bare when they do not.
#[test]
fn counts_csv_quotes_fields_that_need_it() {
    assert_eq!(csv_field("start-A-end"), "start-A-end");
    assert_eq!(csv_field("attB,Tp901"), "\"attB,Tp901\"");
    assert_eq!(csv_field("attB\"Tp901"), "\"attB\"\"Tp901\"");
    assert_eq!(csv_field("a,b\"c"), "\"a,b\"\"c\"");
    // A newline in a name would otherwise end the record early.
    assert_eq!(csv_field("a\nb"), "\"a\nb\"");
    // Quoting reaches the rendered file, not just the helper.
    assert_eq!(
        counts_csv(&["odd,name-A"]),
        "segments,count\n\"odd,name-A\",1\n"
    );
}

// ---------------------------------------------------------------------------
// Equivalence of the linear-space anchor aligner with the library one
// ---------------------------------------------------------------------------

/// The anchor aligner replaced a library implementation, so it must reproduce that
/// aligner's affine gap model exactly. Sequences are drawn from A/T and ACGT because the
/// library aligner compares bases with `==`, which only agrees with the nucleotide-set
/// rule on unambiguous input: with gap_open and gap_extend both -2 a run of k
/// gaps costs -2 - 2k, not -2k. Scores decide whether a read is accepted at all, so they
/// are checked for exact equality over randomised sequences.
///
/// Spans are checked too but not required to match everywhere. On degenerate input -
/// short, low-complexity sequence - several alignments are often exactly co-optimal, and
/// which one gets reported is an arbitrary internal choice that the two implementations
/// make differently. Because the scores always agree, a differing span is never one
/// aligner finding a better alignment than the other, only a different tie broken.
/// `anchor_aligners_agree_on_realistic_reads` covers the case that actually matters.
#[test]
fn anchor_aligner_scores_match_the_library_exactly() {
    let mut rng = Rng(0xA5A5_1234_5678_9ABC);
    let (mut cases, mut score_same, mut span_same) = (0usize, 0usize, 0usize);
    let mut first_score_difference = None;
    for alphabet in [&b"AT"[..], &b"ACGT"[..]] {
        for _ in 0..300 {
            let query_len = 1 + rng.below(18);
            let target_len = 1 + rng.below(50);
            let query = rng.seq(query_len, alphabet);
            let target = rng.seq(target_len, alphabet);
            let theirs = align_segment(&query, &target);
            let ours = best_semiglobal(&query, &target, ANCHOR_SCORING);
            cases += 1;
            if ours.score == theirs.score {
                score_same += 1;
            } else if first_score_difference.is_none() {
                first_score_difference = Some(format!(
                    "query={:?} target={:?}: ours {}, library {}",
                    String::from_utf8_lossy(&query),
                    String::from_utf8_lossy(&target),
                    ours.score,
                    theirs.score
                ));
            }
            if (ours.ystart, ours.yend) == (theirs.ystart, theirs.yend) {
                span_same += 1;
            }
        }
    }
    assert_eq!(
        score_same,
        cases,
        "scores must match exactly; first difference: {}",
        first_score_difference.unwrap_or_default()
    );
    // Loose, only to catch a real regression: co-optimal ties are a small minority even
    // on deliberately degenerate input.
    assert!(
        span_same * 100 / cases >= 80,
        "spans agreed in only {span_same} of {cases} cases, far more divergence than \
         co-optimal ties alone would explain"
    );
}

/// How the two anchor aligners compare on realistic input: a 30 bp anchor sought in
/// reads carrying nanopore-like substitutions and indels. Scores must agree exactly,
/// since they decide whether a read is accepted; spans decide where it is trimmed.
#[test]
fn anchor_aligners_agree_on_realistic_reads() {
    let mut rng = Rng(0x00C0_FFEE_1234_5678);
    let anchor = rng.seq(30, b"ACGT");
    let (mut cases, mut score_same, mut span_same) = (0usize, 0usize, 0usize);
    for _ in 0..400 {
        // Plant a copy of the anchor, damaged the way a nanopore read is damaged.
        let mut planted: Vec<u8> = anchor.clone();
        for _ in 0..rng.below(4) {
            let at = rng.below(planted.len());
            match rng.below(3) {
                0 => planted[at] = b"ACGT"[rng.below(4)],
                1 => {
                    planted.remove(at);
                }
                _ => planted.insert(at, b"ACGT"[rng.below(4)]),
            }
        }
        let read_len = 200 + rng.below(400);
        let mut read = rng.seq(read_len, b"ACGT");
        let at = rng.below(read.len());
        for (offset, base) in planted.iter().enumerate() {
            read.insert(at + offset, *base);
        }

        let theirs = align_segment(&anchor, &read);
        let ours = best_semiglobal(&anchor, &read, ANCHOR_SCORING);
        cases += 1;
        if ours.score == theirs.score {
            score_same += 1;
        }
        if (ours.ystart, ours.yend) == (theirs.ystart, theirs.yend) {
            span_same += 1;
        }
    }
    assert_eq!(score_same, cases, "scores must agree on realistic reads");
    assert_eq!(
        span_same, cases,
        "spans agreed in only {span_same} of {cases} realistic cases"
    );
}

// ---------------------------------------------------------------------------
// Chunked reading
// ---------------------------------------------------------------------------

/// Chunks are bounded by both total bases and read count, and stop short at the end of
/// the file. Whatever the bounds, concatenating the chunks reproduces the whole file in
/// order - which is what lets chunked processing give identical output.
#[test]
fn chunks_respect_their_bounds_and_reassemble_the_whole_file() {
    let dir = TempDir::new().unwrap();
    let path = dir.path().join("reads.fastq");
    let reads: Vec<(String, String)> = (0..25)
        .map(|i| (format!("r{i}"), SEG_A.repeat(1 + i % 3)))
        .collect();
    let borrowed: Vec<(&str, &str)> = reads
        .iter()
        .map(|(n, s)| (n.as_str(), s.as_str()))
        .collect();
    write_fastq(&path, &borrowed);

    for (max_bases, max_reads) in [(20, 100), (1000, 4), (40, 2), (usize::MAX, usize::MAX)] {
        let mut stream = ReadStream::open(&path).unwrap();
        let mut seen = Vec::new();
        let mut chunks = 0;
        loop {
            let chunk = stream.next_chunk(max_bases, max_reads).unwrap();
            if chunk.is_empty() {
                break;
            }
            chunks += 1;
            assert!(
                chunk.len() <= max_reads,
                "chunk of {} exceeds the {max_reads} read bound",
                chunk.len()
            );
            // A chunk may pass the base bound only by taking the read that crosses it.
            let bases: usize = chunk.iter().map(|r| r.seq.len()).sum();
            let without_last: usize = bases - chunk.last().unwrap().seq.len();
            assert!(
                without_last < max_bases,
                "chunk held {without_last} bases before its final read, over the {max_bases} bound"
            );
            seen.extend(chunk);
        }
        assert_eq!(as_pairs(&seen), reads, "bounds {max_bases}/{max_reads}");
        assert!(chunks >= 1);
    }
}

/// Reading in chunks must not change what gets classified. The same reads are processed
/// whole and in small chunks, and both must give the same results in the same order and
/// the same summary counts.
#[test]
fn chunked_processing_matches_processing_everything_at_once() {
    let dir = TempDir::new().unwrap();
    let path = dir.path().join("reads.fastq");
    let reads: Vec<(String, String)> = (0..40)
        .map(|i| {
            let body = match i % 4 {
                0 => format!("{}{}{}{}", JUNK_5, START, SEG_A, END),
                1 => rc(&format!("{}{}{}{}", JUNK_5, START, SEG_B, END)),
                2 => format!("{}{}", SPACER, SEG_A),
                _ => format!("{}{}{}{}{}", START, SEG_A, SPACER, SEG_B, END),
            };
            (format!("r{i}"), body)
        })
        .collect();
    let borrowed: Vec<(&str, &str)> = reads
        .iter()
        .map(|(n, s)| (n.as_str(), s.as_str()))
        .collect();
    write_fastq(&path, &borrowed);

    let segments = seg_map(&anchored_segments());
    let classify_all = |max_bases: usize, max_reads: usize| {
        let mut stream = ReadStream::open(&path).unwrap();
        let mut out = Vec::new();
        let mut summary = RunSummary::default();
        loop {
            let chunk = stream.next_chunk(max_bases, max_reads).unwrap();
            if chunk.is_empty() {
                break;
            }
            let (classified, chunk_summary) =
                process_reads(chunk, &segments, (0, 0), false, true, false).unwrap();
            out.extend(classified.into_iter().map(|r| (r.name, r.segments)));
            summary.merge(chunk_summary);
        }
        (out, summary)
    };

    let (whole, whole_summary) = classify_all(usize::MAX, usize::MAX);
    assert!(
        whole_summary.classified > 0,
        "the fixture should classify something"
    );
    for (max_bases, max_reads) in [(1, 1), (50, 3), (500, 7)] {
        let (chunked, chunked_summary) = classify_all(max_bases, max_reads);
        assert_eq!(
            chunked, whole,
            "results differ at bounds {max_bases}/{max_reads}"
        );
        assert_eq!(
            chunked_summary, whole_summary,
            "summary differs at bounds {max_bases}/{max_reads}"
        );
    }
}

/// Byte progress is measured against the file on disk, so it works for compressed input
/// too, where the decompressed size is not known in advance. The fixture is large enough
/// to span several read buffers, so the count genuinely advances as the file is consumed
/// rather than jumping straight to the end.
#[test]
fn byte_progress_tracks_the_file_for_every_format() {
    let dir = TempDir::new().unwrap();
    // Distinct random sequences, so the gzipped fixture stays far larger than a single
    // read buffer; 3000 copies of one sequence would compress down to a few kilobytes.
    let mut rng = Rng(0x5EED_1234_5678_9ABC);
    let reads: Vec<(String, String)> = (0..3000)
        .map(|i| {
            let seq = String::from_utf8(rng.seq(200, b"ACGT")).unwrap();
            (format!("r{i}"), seq)
        })
        .collect();
    let borrowed: Vec<(&str, &str)> = reads
        .iter()
        .map(|(n, s)| (n.as_str(), s.as_str()))
        .collect();
    let bam_reads: Vec<(&str, &str, Flags)> = borrowed
        .iter()
        .map(|(n, s)| (*n, *s, Flags::UNMAPPED))
        .collect();

    type WriteFixture<'a> = Box<dyn Fn(&Path) + 'a>;
    let writers: [(&str, WriteFixture); 3] = [
        (
            "reads.fastq",
            Box::new(|p: &Path| write_fastq(p, &borrowed)),
        ),
        (
            "reads.fastq.gz",
            Box::new(|p: &Path| write_fastq_gz(p, &borrowed)),
        ),
        ("reads.bam", Box::new(|p: &Path| write_bam(p, &bam_reads))),
    ];
    for (name, write_fixture) in writers {
        let path = dir.path().join(name);
        write_fixture(&path);
        let size = std::fs::metadata(&path).unwrap().len();

        let mut stream = ReadStream::open(&path).unwrap();
        let mut positions = vec![stream.bytes_read()];
        let mut total_reads = 0;
        loop {
            let chunk = stream.next_chunk(2000, 200).unwrap();
            if chunk.is_empty() {
                break;
            }
            total_reads += chunk.len();
            positions.push(stream.bytes_read());
        }

        assert_eq!(total_reads, reads.len(), "{name}: not every read came back");
        assert!(
            positions.windows(2).all(|w| w[1] >= w[0]),
            "{name}: progress went backwards: {positions:?}"
        );
        let consumed = *positions.last().unwrap();
        assert!(
            consumed > positions[0],
            "{name}: reading advanced the count no further than opening did"
        );
        assert!(
            consumed <= size,
            "{name}: counted {consumed} bytes from a {size} byte file"
        );
        // Progress should arrive in steps, not a single jump at the very end.
        assert!(
            positions.len() > 3,
            "{name}: expected several chunks, got {}",
            positions.len() - 1
        );
    }
}

// ---------------------------------------------------------------------------
// Detailed output
// ---------------------------------------------------------------------------

/// Pull the spans out of a located string. Entries are comma separated and the two
/// positions within one are separated by a colon, so a plain split on ',' is enough -
/// which is the point of the colon, and is asserted on directly below.
fn located_spans(located: &str) -> Vec<(String, usize, usize)> {
    if located.is_empty() {
        return Vec::new();
    }
    located
        .split(',')
        .map(|entry| {
            let (name, span) = entry
                .split_once('[')
                .expect("an entry looks like name[a:b]");
            let (start, end) = span
                .trim_end_matches(']')
                .split_once(':')
                .expect("positions are separated by a colon");
            (
                name.to_string(),
                start.parse().unwrap(),
                end.parse().unwrap(),
            )
        })
        .collect()
}

/// The separator between entries and the separator within one are different characters,
/// so splitting the column on ',' can never cut an entry in half.
#[test]
fn entries_are_separated_differently_from_the_positions_within_them() {
    let read = format!("{}{}{}{}", SEG_A, SPACER, SEG_B, SEG_A);
    let (_, located, _) = run_detailed(&read, &[("A", SEG_A), ("B", SEG_B)], false, false);
    assert_eq!(located, "A[1:20],B[33:52],A[53:72]");
    assert_eq!(
        located.split(',').count(),
        3,
        "a plain comma split yields one piece per segment: {located}"
    );
}

/// Without the flag there is no detail to write, so the extra columns cost nothing on a
/// run that does not ask for them.
#[test]
fn the_extra_columns_are_absent_unless_asked_for() {
    let read = format!("{}{}{}", SPACER, SEG_A, SPACER);
    let (classified, _) = process_reads(
        vec![raw("read", &read)],
        &seg_map(&[("A", SEG_A)]),
        (0, 0),
        false,
        false,
        false,
    )
    .unwrap();
    assert!(classified[0].detail.is_none());
    assert_eq!(classified[0].segments, "A");
}

/// The positions index the sequence written beside them, counting from 1 and including
/// both ends. Cutting the sequence at the span has to give back the segment.
#[test]
fn positions_are_one_based_and_index_the_sequence_column() {
    let read = format!("{}{}{}", SPACER, SEG_A, SPACER);
    let (segments, located, sequence) = run_detailed(&read, &[("A", SEG_A)], false, false);

    assert_eq!(segments, "A");
    // SPACER is 12 bp, so a 20 bp segment planted after it runs from 13 to 32.
    assert_eq!(located, "A[13:32]");
    let (_, start, end) = located_spans(&located)[0].clone();
    assert_eq!(&sequence[start - 1..end], SEG_A, "the span cuts out SEG_A");
    assert_eq!(end - start + 1, SEG_A.len(), "inclusive of both ends");
}

/// One entry per occurrence, in read order, comma separated - matching the segment
/// string beside it name for name.
#[test]
fn every_occurrence_gets_its_own_entry_in_read_order() {
    let read = format!("{}{}{}{}", SEG_A, SPACER, SEG_B, SEG_A);
    let (segments, located, sequence) =
        run_detailed(&read, &[("A", SEG_A), ("B", SEG_B)], false, false);

    assert_eq!(segments, "A-B-A");
    assert_eq!(located, "A[1:20],B[33:52],A[53:72]");
    let spans = located_spans(&located);
    assert_eq!(
        spans
            .iter()
            .map(|(name, ..)| name.as_str())
            .collect::<Vec<_>>(),
        segments.split('-').collect::<Vec<_>>(),
        "the two columns name the same segments in the same order"
    );
    for (name, start, end) in spans {
        let expected = if name == "A" { SEG_A } else { SEG_B };
        assert_eq!(&sequence[start - 1..end], expected, "span for {name}");
    }
}

/// A segment found on the reverse strand carries the same trailing '*' it has in the
/// segment string, and its span still reads along the extracted sequence.
#[test]
fn reverse_strand_segments_keep_their_star() {
    let read = format!("{}{}{}", SPACER, rc(SEG_A), SPACER);
    let (segments, located, sequence) = run_detailed(&read, &[("A", SEG_A)], false, false);

    assert_eq!(segments, "A*");
    assert_eq!(located, "A*[13:32]");
    assert_eq!(
        &sequence[12..32],
        rc(SEG_A),
        "the span cuts out the revcomp"
    );
}

/// A read with nothing above threshold has an empty located column rather than a missing
/// one, so every row still has the same number of columns.
#[test]
fn a_read_with_no_segments_has_an_empty_located_column() {
    let read = SPACER.repeat(3);
    let (segments, located, sequence) = run_detailed(&read, &[("A", SEG_A)], false, false);
    assert_eq!((segments.as_str(), located.as_str()), ("", ""));
    assert_eq!(sequence, read, "the sequence is still written");
}

/// Without anchors the extracted sequence is the read exactly as it arrived, so the
/// positions are positions in the read itself.
#[test]
fn the_sequence_is_the_whole_read_when_there_are_no_anchors() {
    let read = format!("{}{}{}", JUNK_5, SEG_A, JUNK_3);
    let (_, located, sequence) = run_detailed(&read, &[("A", SEG_A)], false, false);
    assert_eq!(sequence, read);
    assert_eq!(
        located, "A[11:30]",
        "JUNK_5 is 10 bp, so SEG_A starts at 11"
    );
}

/// With anchors the extracted sequence is the trimmed span between them, and the
/// positions count along that rather than along the read as it arrived.
#[test]
fn the_sequence_is_trimmed_to_the_anchors_and_positions_follow_it() {
    let read = format!("{}{}{}", JUNK_5, construct(), JUNK_3);
    let (segments, located, sequence) = run_detailed(&read, &anchored_segments(), false, true);

    assert_eq!(segments, "start-A-B-end");
    assert_eq!(
        sequence,
        construct(),
        "trimmed to the anchors, junk removed"
    );
    // Positions along the construct, not the read: START is 20 bp, so A begins at 21.
    assert_eq!(located, "start[1:20],A[21:40],B[53:72],end[73:92]");
    assert_eq!(&sequence[20..40], SEG_A);
    assert_eq!(&sequence[52..72], SEG_B);
}

/// The anchors are located like everything else, so the two columns match name for name
/// and the anchors bound the sequence exactly: 'start' opens it and 'end' closes it.
/// Reporting them is what lets a run be diagnosed from the output alone - a trimmed read
/// that looks wrong can be traced to where the anchors were actually found.
#[test]
fn the_anchors_are_located_like_any_other_segment() {
    let read = format!("{}{}{}", JUNK_5, construct(), JUNK_3);
    let (segments, located, sequence) = run_detailed(&read, &anchored_segments(), false, true);

    let spans = located_spans(&located);
    assert_eq!(
        spans
            .iter()
            .map(|(name, ..)| name.as_str())
            .collect::<Vec<_>>(),
        segments.split('-').collect::<Vec<_>>(),
        "every name in the segment string has a span beside it"
    );
    let (first, last) = (spans.first().unwrap(), spans.last().unwrap());
    assert_eq!((first.0.as_str(), first.1), ("start", 1), "start opens it");
    assert_eq!(
        (last.0.as_str(), last.2),
        ("end", sequence.len()),
        "end closes it"
    );
    // ...and the anchor spans cut the anchors back out.
    assert_eq!(&sequence[first.1 - 1..first.2], START);
    assert_eq!(&sequence[last.1 - 1..last.2], END);
}

/// An anchor that aligned across a different number of bases than it is long reports the
/// span it actually covered, not its nominal length. That is the point of reporting them:
/// a degraded anchor shows up here rather than having to be inferred.
#[test]
fn an_anchor_span_reflects_the_alignment_not_the_anchor_length() {
    // A read whose start anchor is missing its first base, flush with the read start so
    // there is nothing before it for the aligner to absorb that base against.
    let read = format!("{}{}{}", &START[1..], SEG_A, END);
    let (segments, located, sequence) = run_detailed(&read, &anchored_segments(), false, true);

    assert_eq!(segments, "start-A-end");
    let spans = located_spans(&located);
    assert_eq!(spans[0], ("start".to_string(), 1, 19), "19 bases, not 20");
    assert_eq!(sequence.len(), 19 + SEG_A.len() + END.len());
    assert_eq!(&sequence[..19], &START[1..]);
}

/// A read sequenced the other way round is reoriented before classification, so the
/// sequence column comes out forward and the positions match the forward read exactly.
/// This is the case where the positions would be wrong if they were taken from the read
/// as it arrived rather than from the sequence written beside them.
#[test]
fn a_reverse_read_is_reoriented_before_its_positions_are_taken() {
    let forward = format!("{}{}{}", JUNK_5, construct(), JUNK_3);
    let reverse = rc(&forward);
    assert_eq!(
        run_detailed(&reverse, &anchored_segments(), false, true),
        run_detailed(&forward, &anchored_segments(), false, true),
        "both strands describe the same construct"
    );
}

/// A read that wraps the origin is rotated back into order first, so with --circular the
/// positions describe the rotated sequence that is written out.
#[test]
fn a_wrapped_read_is_rotated_before_its_positions_are_taken() {
    let (_, located, sequence) =
        run_detailed(&rotated_construct(), &anchored_segments(), true, true);
    assert_eq!(sequence, construct(), "rotated back to start-...-end");
    assert_eq!(located, "start[1:20],A[21:40],B[53:72],end[73:92]");
}

// ---------------------------------------------------------------------------
// Report formatting
// ---------------------------------------------------------------------------

/// The calendar arithmetic behind the timestamps, checked against dates whose answers
/// are known independently. Leap years and century rules are where a hand-rolled
/// conversion goes wrong, so those are what is pinned.
#[test]
fn days_since_the_epoch_convert_to_the_right_calendar_date() {
    for (days, expected) in [
        (0, (1970, 1, 1)),
        (1, (1970, 1, 2)),
        (30, (1970, 1, 31)),
        (31, (1970, 2, 1)),
        (364, (1970, 12, 31)),
        (365, (1971, 1, 1)),
        // 1972 was a leap year, so it has a 29 February.
        (789, (1972, 2, 29)),
        (790, (1972, 3, 1)),
        // 2000 was a leap year despite being a century: divisible by 400.
        (11016, (2000, 2, 29)),
        // 1900 was not, being a century that is not divisible by 400: 28 February is
        // followed straight by 1 March.
        (-25509, (1900, 2, 28)),
        (-25508, (1900, 3, 1)),
        (20000, (2024, 10, 4)),
    ] {
        assert_eq!(
            civil_from_days(days),
            expected,
            "{days} days since the epoch"
        );
    }
}

/// Timestamps come out in a fixed, sortable form, and a time the clock cannot express
/// says so rather than panicking.
#[test]
fn timestamps_are_formatted_as_sortable_utc() {
    let at = |secs: u64| format_utc(SystemTime::UNIX_EPOCH + Duration::from_secs(secs));
    assert_eq!(at(0), "1970-01-01 00:00:00 UTC");
    assert_eq!(at(86_399), "1970-01-01 23:59:59 UTC");
    assert_eq!(at(86_400), "1970-01-02 00:00:00 UTC");
    assert_eq!(at(1_700_000_000), "2023-11-14 22:13:20 UTC");
    // Sub-second precision is dropped rather than rounded up into the next second.
    assert_eq!(
        format_utc(SystemTime::UNIX_EPOCH + Duration::from_millis(999)),
        "1970-01-01 00:00:00 UTC"
    );
    // A clock set before 1970 is not worth handling, but must not bring the run down
    // after the work is already done.
    assert_eq!(
        format_utc(SystemTime::UNIX_EPOCH - Duration::from_secs(1)),
        "unknown"
    );
}

/// Sizes are quoted at a scale a person can read, with whole bytes left whole.
#[test]
fn byte_counts_are_scaled_for_reading() {
    assert_eq!(human_bytes(0), "0 B");
    assert_eq!(human_bytes(1023), "1023 B");
    assert_eq!(human_bytes(1024), "1.0 KiB");
    assert_eq!(human_bytes(1536), "1.5 KiB");
    assert_eq!(human_bytes(1024 * 1024), "1.0 MiB");
    assert_eq!(human_bytes(3 * 1024 * 1024 * 1024), "3.0 GiB");
    // Beyond the largest unit it keeps scaling in that unit rather than wrapping round.
    assert!(human_bytes(u64::MAX).ends_with(" TiB"));
}

/// Base counts scale by powers of ten, which is how sequence lengths are always quoted -
/// unlike file sizes, which go in powers of two.
#[test]
fn base_counts_are_scaled_by_powers_of_ten() {
    assert_eq!(human_bases(0.0), "0.0 base");
    assert_eq!(human_bases(999.0), "999.0 base");
    assert_eq!(human_bases(1_000.0), "1.0 kbase");
    assert_eq!(human_bases(7_561_380.0), "7.6 Mbase");
    assert_eq!(human_bases(2.5e9), "2.5 Gbase");
}

/// Durations stay in seconds while that reads well, then break into minutes and hours.
#[test]
fn durations_are_scaled_for_reading() {
    let secs = |s: f64| human_duration(Duration::from_secs_f64(s));
    assert_eq!(secs(0.0), "0.00 s");
    assert_eq!(secs(0.314), "0.31 s");
    assert_eq!(secs(59.99), "59.99 s");
    assert_eq!(secs(60.0), "1 m 00 s");
    assert_eq!(secs(3599.0), "59 m 59 s");
    assert_eq!(secs(3600.0), "1 h 00 m 00 s");
    assert_eq!(secs(7384.0), "2 h 03 m 04 s");
}

/// A run too fast for the clock to separate has no meaningful rate, so the rates are
/// left out rather than printed as an infinity.
#[test]
fn a_run_with_no_elapsed_time_quotes_no_rate() {
    let report = RunReport {
        started: SystemTime::UNIX_EPOCH,
        finished: SystemTime::UNIX_EPOCH,
        elapsed: Duration::ZERO,
        segments_file: "refs.fasta".to_string(),
        segments_loaded: 2,
        sequences_file: "reads.fastq".to_string(),
        sequences_format: "FASTQ",
        sequences_bytes: 0,
        reads: 0,
        bases: 0,
        options: Vec::new(),
    }
    .render(&RunSummary::default());
    assert!(!report.contains("per second"), "{report}");
    assert!(
        !report.contains("inf") && !report.contains("NaN"),
        "{report}"
    );
    // Everything that does not depend on the clock is still there.
    assert!(
        report.contains("refs.fasta") && report.contains("FASTQ"),
        "{report}"
    );
}

/// The report is one document: the sections share a label column, so it reads as a
/// table rather than as several tables stacked up.
#[test]
fn the_report_sections_share_one_label_column() {
    let report = RunReport {
        started: SystemTime::UNIX_EPOCH,
        finished: SystemTime::UNIX_EPOCH + Duration::from_secs(2),
        elapsed: Duration::from_secs(2),
        segments_file: "refs.fasta".to_string(),
        segments_loaded: 2,
        sequences_file: "reads.fastq".to_string(),
        sequences_format: "FASTQ",
        sequences_bytes: 4096,
        reads: 10,
        bases: 1000,
        options: vec![("--threads".to_string(), "4".to_string())],
    }
    .render(&RunSummary {
        classified: 8,
        too_short: 2,
        ..Default::default()
    });

    // Values start where counts end: one column for the whole document.
    let starts: Vec<usize> = report
        .lines()
        .filter(|line| line.starts_with("  ") && line.contains("  ") && !line.trim().is_empty())
        .filter_map(|line| line.rfind("  ").map(|gap| gap + 2))
        .collect();
    assert!(
        starts.len() > 8,
        "expected rows from every section:\n{report}"
    );
    // Each section is introduced by a heading on its own line.
    for heading in ["Run", "Input", "Options", "Summary", "Throughput"] {
        assert!(
            report.lines().any(|line| line == heading),
            "missing the {heading} heading:\n{report}"
        );
    }
}

// ---------------------------------------------------------------------------
// Summary layout
// ---------------------------------------------------------------------------

/// A summary line with its trailing percentage removed, so what is left ends with the
/// count. Used by the alignment tests below to find where each column sits.
fn without_percentage(line: &str) -> &str {
    match line.rfind('(') {
        Some(open) if line.ends_with("%)") => line[..open].trim_end(),
        _ => line.trim_end(),
    }
}

/// The summary lines that carry a count, i.e. everything but the heading.
fn counted_lines(report: &str) -> Vec<&str> {
    report
        .lines()
        .filter(|line| line.starts_with(' '))
        .collect()
}

/// Every count ends in the same column whatever its nesting depth. The indented
/// breakdown rows used to be laid out to their own hardcoded width, which put them
/// nineteen columns to the right of the totals they were breaking down.
#[test]
fn every_count_ends_in_the_same_column() {
    let report = RunSummary {
        classified: 8,
        too_short: 2,
        anchors_out_of_order: 250,
        unreadable: 9,
        ..Default::default()
    }
    .render();
    let columns: Vec<usize> = counted_lines(&report)
        .iter()
        .map(|line| without_percentage(line).len())
        .collect();
    assert_eq!(columns.len(), 6, "expected every row to count:\n{report}");
    assert!(
        columns.iter().all(|column| *column == columns[0]),
        "counts do not share a column:\n{report}"
    );
}

/// Percentages are right-aligned as well, so the decimal point and the closing bracket
/// line up whether the figure is 0.0% or 100.0% - two characters wider.
#[test]
fn every_percentage_ends_in_the_same_column() {
    let report = RunSummary {
        classified: 20,
        ..Default::default()
    }
    .render();
    assert!(report.contains("(100.0%)"), "{report}");
    assert!(report.contains("(0.0%)"), "{report}");
    let ends: Vec<usize> = report
        .lines()
        .filter(|line| line.ends_with("%)"))
        .map(str::len)
        .collect();
    assert_eq!(ends.len(), 2, "both shares are quoted:\n{report}");
    assert_eq!(
        ends[0], ends[1],
        "percentages do not share a column:\n{report}"
    );
}

/// The columns are measured from the rows being printed rather than fixed, so a small
/// run is not padded out to the width of a large one, and a large one is not squeezed.
#[test]
fn columns_are_sized_to_the_numbers_they_hold() {
    let column_of = |classified: usize| {
        let report = RunSummary {
            classified,
            ..Default::default()
        }
        .render();
        without_percentage(counted_lines(&report)[0]).len()
    };
    assert!(
        column_of(7_000_000) > column_of(7),
        "the column should grow with the numbers in it"
    );
    let large = RunSummary {
        classified: 7_000_000,
        ..Default::default()
    }
    .render();
    assert!(large.contains("7000000"), "nothing is truncated:\n{large}");
}

/// However long the label, there is always a gap before its count. The widest row is the
/// one that sets the column, and it is the one at risk of the two running together.
#[test]
fn no_label_runs_into_its_count() {
    let report = RunSummary {
        classified: 1,
        anchors_on_opposite_strands: 1,
        unreadable: 1,
        ..Default::default()
    }
    .render();
    for line in counted_lines(&report) {
        let upto_count = without_percentage(line);
        let count_start = upto_count.rfind(' ').map_or(0, |space| space + 1);
        assert!(
            upto_count[..count_start].ends_with("  "),
            "label runs into its count: {line:?}"
        );
    }
}

/// With no reads there is no share to quote, so the percentage column disappears rather
/// than dividing by zero, and what is left still lines up.
#[test]
fn a_summary_of_no_reads_has_no_percentage_column() {
    let report = RunSummary::default().render();
    assert!(!report.contains('%'), "no percentages to show:\n{report}");
    assert!(!report.contains("NaN"), "{report}");
    let columns: Vec<usize> = counted_lines(&report).iter().map(|l| l.len()).collect();
    assert!(
        columns.iter().all(|column| *column == columns[0]),
        "counts do not share a column:\n{report}"
    );
}
