// Segment alignment tool
// Author: Thomas E. Gorochowski <tom@chofski.co.uk>

use bio::alphabets::dna;
use bio::io::fastq::FastqRead;
use bio::io::{fasta, fastq};
use clap::Parser;
use indicatif::{ProgressBar, ProgressStyle};
use noodles_bam as bam;
use noodles_bgzf as bgzf;
use rayon::prelude::*;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Read, Write};
use std::path::Path;

#[cfg(test)]
mod tests;

/// Every fallible routine reports a plain message meant to be read by the user, which
/// `main` prints and exits on. Nothing in normal operation should reach a panic.
type Result<T> = std::result::Result<T, String>;

/// Report a recoverable problem without stopping the run.
fn warn(message: String) {
    eprintln!("warning: {message}");
}

/// Read at most this much sequence into memory at once. This bounds peak memory
/// independently of how large the input file is, while staying big enough to keep every
/// thread busy and to update the progress bar often.
const MAX_CHUNK_BASES: usize = 8 << 20;
/// ...and at most this many reads, so a file of very short reads still chunks sensibly.
const MAX_CHUNK_READS: usize = 100_000;

/// Give up on a file after this many malformed records in a row. A file that never
/// yields another usable record is a broken file, not a file with a few bad reads.
const MAX_CONSECUTIVE_BAD_RECORDS: usize = 100;

/// Loads a FASTA file containing the `start` and `end` sequences plus and
/// other segment sequences that should be used when classifying a read.
///
/// Sequences are upper-cased on the way in so that a soft-masked or lower-case FASTA
/// still matches reads. Every problem here is fatal: a segment that cannot be used is
/// a mistake in the input, and carrying on would silently change what gets reported.
fn load_fasta(filename: &Path) -> Result<HashMap<String, Vec<u8>>> {
    let mut fasta_data: HashMap<String, Vec<u8>> = HashMap::new();
    // Opened directly rather than via `from_file` so the message carries the actual
    // reason (e.g. "No such file or directory") instead of a generic wrapper.
    let file = File::open(filename)
        .map_err(|e| format!("Could not open segments file '{}': {e}", filename.display()))?;
    let reader = fasta::Reader::new(file);
    for (index, result) in reader.records().enumerate() {
        let record = result.map_err(|e| {
            format!(
                "Could not parse segment {} in '{}': {e}",
                index + 1,
                filename.display()
            )
        })?;
        let seg_name = record.id().to_string();
        let seg_seq = record.seq().to_ascii_uppercase();
        if seg_seq.is_empty() {
            return Err(format!(
                "Segment '{seg_name}' in '{}' has no sequence. Every segment needs at least one base.",
                filename.display()
            ));
        }
        if fasta_data.contains_key(&seg_name) {
            return Err(format!(
                "Segment '{seg_name}' is defined more than once in '{}'. Segment names must be unique, \
                 otherwise it is ambiguous which sequence should be searched for.",
                filename.display()
            ));
        }
        fasta_data.insert(seg_name, seg_seq);
    }
    if fasta_data.is_empty() {
        return Err(format!(
            "No segments found in '{}'. Expected a FASTA file with one record per segment.",
            filename.display()
        ));
    }
    Ok(fasta_data)
}

/// Structure to hold the key sequence and segment information about a read.
struct SegmentedRead {
    name: String,
    segments: String,
}

/// A read name/sequence pair, independent of the input file format it came from.
struct RawRead {
    name: String,
    seq: Vec<u8>,
}

/// Reads loaded from a sequences file, and how many records had to be skipped to get
/// them (malformed, unnamed, or carrying no sequence).
#[cfg(test)]
struct LoadedReads {
    reads: Vec<RawRead>,
    #[allow(dead_code)] // mirrors what `run` reports; kept so the shape matches
    skipped: usize,
}

/// Why a read produced no classification. Every read that is not classified falls into
/// exactly one of these, so the counts in the summary always add up.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Rejected {
    TooShort,
    TooLong,
    StartAnchorNotFound,
    EndAnchorNotFound,
    NeitherAnchorFound,
    AnchorsOnOppositeStrands,
    AnchorsOutOfOrder,
}

/// The result of trying to classify one read.
enum Outcome {
    Classified(SegmentedRead),
    Rejected(Rejected),
}

/// Tally of what happened to every read, reported at the end of the run.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct RunSummary {
    /// Records skipped while reading the input file, before classification.
    pub unreadable: usize,
    pub classified: usize,
    pub too_short: usize,
    pub too_long: usize,
    pub start_anchor_not_found: usize,
    pub end_anchor_not_found: usize,
    pub neither_anchor_found: usize,
    pub anchors_on_opposite_strands: usize,
    pub anchors_out_of_order: usize,
}

impl RunSummary {
    /// Record one rejected read under its reason.
    fn count(&mut self, reason: Rejected) {
        let slot = match reason {
            Rejected::TooShort => &mut self.too_short,
            Rejected::TooLong => &mut self.too_long,
            Rejected::StartAnchorNotFound => &mut self.start_anchor_not_found,
            Rejected::EndAnchorNotFound => &mut self.end_anchor_not_found,
            Rejected::NeitherAnchorFound => &mut self.neither_anchor_found,
            Rejected::AnchorsOnOppositeStrands => &mut self.anchors_on_opposite_strands,
            Rejected::AnchorsOutOfOrder => &mut self.anchors_out_of_order,
        };
        *slot += 1;
    }

    /// Fold another chunk's tally into this one.
    fn merge(&mut self, other: RunSummary) {
        self.unreadable += other.unreadable;
        self.classified += other.classified;
        self.too_short += other.too_short;
        self.too_long += other.too_long;
        self.start_anchor_not_found += other.start_anchor_not_found;
        self.end_anchor_not_found += other.end_anchor_not_found;
        self.neither_anchor_found += other.neither_anchor_found;
        self.anchors_on_opposite_strands += other.anchors_on_opposite_strands;
        self.anchors_out_of_order += other.anchors_out_of_order;
    }

    /// Reads that reached classification but produced no result.
    pub fn not_classified(&self) -> usize {
        self.too_short
            + self.too_long
            + self.start_anchor_not_found
            + self.end_anchor_not_found
            + self.neither_anchor_found
            + self.anchors_on_opposite_strands
            + self.anchors_out_of_order
    }

    /// Reads that reached classification, whatever the outcome.
    pub fn reads(&self) -> usize {
        self.classified + self.not_classified()
    }

    /// Render the end-of-run summary. Categories that cannot arise for the given
    /// settings - the anchor ones without --start-end-segs, for instance - stay at zero
    /// and are left out, so the report only shows what actually happened.
    fn render(&self) -> String {
        let pct = |n: usize| {
            if self.reads() == 0 {
                String::new()
            } else {
                format!("  ({:.1}%)", 100.0 * n as f64 / self.reads() as f64)
            }
        };
        let mut out = String::from("\nSummary\n");
        out.push_str(&format!("  reads read:      {:>9}\n", self.reads()));
        out.push_str(&format!(
            "  classified:      {:>9}{}\n",
            self.classified,
            pct(self.classified)
        ));
        out.push_str(&format!(
            "  not classified:  {:>9}{}\n",
            self.not_classified(),
            pct(self.not_classified())
        ));
        for (label, count) in [
            ("read too short", self.too_short),
            ("read too long", self.too_long),
            ("'start' segment not found", self.start_anchor_not_found),
            ("'end' segment not found", self.end_anchor_not_found),
            ("neither segment found", self.neither_anchor_found),
            (
                "'start'/'end' on opposite strands",
                self.anchors_on_opposite_strands,
            ),
            ("'start'/'end' out of order", self.anchors_out_of_order),
        ] {
            if count > 0 {
                out.push_str(&format!("    {label:<36}{count:>7}\n"));
            }
        }
        if self.unreadable > 0 {
            out.push_str(&format!(
                "  records skipped while reading the input: {}\n",
                self.unreadable
            ));
        }
        out
    }
}

/// The auto-detected format of a sequences file.
enum InputFormat {
    Fastq,
    FastqGz,
    Bam,
}

/// Inspect a file to determine whether it's FASTQ (starts with '@'), gzipped FASTQ, or BAM.
/// BAM is always BGZF-compressed (i.e. gzip-wrapped), so a gzip magic number alone doesn't
/// distinguish it from a gzipped FASTQ file - the decompressed content is peeked to tell
/// them apart (BAM starts with the "BAM\1" magic bytes, FASTQ starts with '@').
fn detect_format(filename: &Path) -> Result<InputFormat> {
    let open = |what: &str| {
        File::open(filename).map_err(|e| {
            format!(
                "Could not open sequences file '{}'{what}: {e}",
                filename.display()
            )
        })
    };
    let mut f = open("")?;
    let mut magic = [0u8; 2];
    let n = f.read(&mut magic).map_err(|e| {
        format!(
            "Could not read sequences file '{}': {e}",
            filename.display()
        )
    })?;
    if n == 2 && magic[0] == 0x1f && magic[1] == 0x8b {
        let mut decoder = flate2::read::MultiGzDecoder::new(open(" for decompression")?);
        let mut peek = [0u8; 4];
        let peek_n = decoder.read(&mut peek).map_err(|e| {
            format!(
                "Could not decompress sequences file '{}': {e}. It starts like a gzip file but \
                 does not decompress; it may be truncated or corrupt.",
                filename.display()
            )
        })?;
        if peek_n == 4 && &peek == b"BAM\x01" {
            Ok(InputFormat::Bam)
        } else if peek_n >= 1 && peek[0] == b'@' {
            Ok(InputFormat::FastqGz)
        } else {
            Err(format!(
                "Could not determine the format of '{}'. It is gzip-compressed but its contents are \
                 neither FASTQ (records start with '@') nor BAM.",
                filename.display()
            ))
        }
    } else if n >= 1 && magic[0] == b'@' {
        Ok(InputFormat::Fastq)
    } else {
        Err(format!(
            "Could not determine the format of '{}'. Expected FASTQ (records start with '@'), \
             gzipped FASTQ, or an unaligned BAM.",
            filename.display()
        ))
    }
}

/// Counts bytes pulled from the file underneath, so progress can be reported against the
/// size on disk. Because the count sits below any decompression, the same measure works
/// for plain, gzipped and BAM input alike: how much of the file has been consumed.
struct CountingReader<R> {
    inner: R,
    bytes: std::sync::Arc<std::sync::atomic::AtomicU64>,
}

impl<R: Read> Read for CountingReader<R> {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        let n = self.inner.read(buf)?;
        self.bytes
            .fetch_add(n as u64, std::sync::atomic::Ordering::Relaxed);
        Ok(n)
    }
}

/// The parser behind a `ReadStream`, one per accepted input format.
enum Source {
    Fastq(fastq::Reader<std::io::BufReader<Box<dyn Read>>>),
    Bam(Box<bam::io::Reader<bgzf::io::Reader<CountingReader<File>>>>),
}

/// Streams reads from a sequences file in bounded chunks, so peak memory follows the
/// chunk size rather than the size of the input.
///
/// Unusable records are skipped here rather than aborting the run, and the tally is
/// reported once the whole file has been read.
struct ReadStream {
    source: Source,
    filename: std::path::PathBuf,
    bytes_read: std::sync::Arc<std::sync::atomic::AtomicU64>,
    record: fastq::Record,
    bam_record: bam::Record,
    skipped: usize,
    index: usize,
    consecutive_errors: usize,
}

impl ReadStream {
    /// Open a sequences file, detecting its format from the contents.
    fn open(filename: &Path) -> Result<Self> {
        let bytes_read = std::sync::Arc::new(std::sync::atomic::AtomicU64::new(0));
        let counted =
            |bytes: &std::sync::Arc<std::sync::atomic::AtomicU64>| -> Result<CountingReader<File>> {
                let file = File::open(filename).map_err(|e| {
                    format!(
                        "Could not open sequences file '{}': {e}",
                        filename.display()
                    )
                })?;
                Ok(CountingReader {
                    inner: file,
                    bytes: bytes.clone(),
                })
            };
        let source = match detect_format(filename)? {
            InputFormat::Fastq => {
                let raw: Box<dyn Read> = Box::new(counted(&bytes_read)?);
                Source::Fastq(fastq::Reader::new(raw))
            }
            InputFormat::FastqGz => {
                let raw: Box<dyn Read> =
                    Box::new(flate2::read::MultiGzDecoder::new(counted(&bytes_read)?));
                Source::Fastq(fastq::Reader::new(raw))
            }
            InputFormat::Bam => {
                let mut reader = bam::io::Reader::new(counted(&bytes_read)?);
                reader.read_header().map_err(|e| {
                    format!(
                        "Could not read the BAM header of '{}': {e}",
                        filename.display()
                    )
                })?;
                Source::Bam(Box::new(reader))
            }
        };
        Ok(ReadStream {
            source,
            filename: filename.to_path_buf(),
            bytes_read,
            record: fastq::Record::new(),
            bam_record: bam::Record::default(),
            skipped: 0,
            index: 0,
            consecutive_errors: 0,
        })
    }

    /// Bytes of the file consumed so far.
    fn bytes_read(&self) -> u64 {
        self.bytes_read.load(std::sync::atomic::Ordering::Relaxed)
    }

    /// Take up to `max_reads` reads totalling at most `max_bases`, or fewer at the end
    /// of the file. An empty chunk means there is nothing left.
    fn next_chunk(&mut self, max_bases: usize, max_reads: usize) -> Result<Vec<RawRead>> {
        let mut chunk = Vec::new();
        let mut bases = 0usize;
        while chunk.len() < max_reads && bases < max_bases {
            match self.next_read()? {
                Some(read) => {
                    bases += read.seq.len();
                    chunk.push(read);
                }
                None => break,
            }
        }
        Ok(chunk)
    }

    /// The next usable read, skipping over any that cannot be used.
    fn next_read(&mut self) -> Result<Option<RawRead>> {
        match &mut self.source {
            Source::Fastq(reader) => loop {
                match reader.read(&mut self.record) {
                    Ok(()) => {
                        if self.record.is_empty() {
                            return Ok(None);
                        }
                        self.index += 1;
                        self.consecutive_errors = 0;
                        if self.record.seq().is_empty() {
                            warn(format!(
                                "skipping read '{}' in '{}': it has no sequence",
                                self.record.id(),
                                self.filename.display()
                            ));
                            self.skipped += 1;
                            continue;
                        }
                        return Ok(Some(RawRead {
                            name: self.record.id().to_string(),
                            seq: self.record.seq().to_ascii_uppercase(),
                        }));
                    }
                    Err(e) => {
                        self.index += 1;
                        self.skipped += 1;
                        self.consecutive_errors += 1;
                        // The name is only known if the parser reached the header line.
                        if self.record.id().is_empty() {
                            warn(format!(
                                "skipping malformed record {} in '{}': {e}",
                                self.index,
                                self.filename.display()
                            ));
                        } else {
                            warn(format!(
                                "skipping malformed read '{}' in '{}': {e}",
                                self.record.id(),
                                self.filename.display()
                            ));
                        }
                        if self.consecutive_errors >= MAX_CONSECUTIVE_BAD_RECORDS {
                            return Err(format!(
                                "Gave up on '{}' after {} malformed records in a row. The file \
                                 is probably truncated or not FASTQ at all.",
                                self.filename.display(),
                                self.consecutive_errors
                            ));
                        }
                    }
                }
            },
            Source::Bam(reader) => loop {
                // A record that fails to decode leaves the position in the compressed
                // stream unknown, so unlike FASTQ there is no safe way to continue.
                let read_bytes = reader.read_record(&mut self.bam_record).map_err(|e| {
                    format!(
                        "Could not read BAM record {} in '{}': {e}. The file may be truncated \
                         or corrupt.",
                        self.index + 1,
                        self.filename.display()
                    )
                })?;
                if read_bytes == 0 {
                    return Ok(None);
                }
                self.index += 1;
                let flags = self.bam_record.flags();
                let name = self.bam_record.name().map(|n| n.to_string());
                if !flags.is_unmapped() || flags.is_secondary() || flags.is_supplementary() {
                    return Err(format!(
                        "'{}' is an aligned BAM file: record {} ('{}') is mapped to a reference \
                         or is a secondary/supplementary alignment. segment needs an unaligned \
                         BAM, because an aligned one repeats reads across alignments and \
                         hard-clips their sequence. Convert it first, for example with \
                         `samtools fastq in.bam > out.fastq`.",
                        self.filename.display(),
                        self.index,
                        name.as_deref().unwrap_or("unnamed"),
                    ));
                }
                let Some(name) = name else {
                    warn(format!(
                        "skipping BAM record {} in '{}': it has no read name",
                        self.index,
                        self.filename.display()
                    ));
                    self.skipped += 1;
                    continue;
                };
                let seq: Vec<u8> = self.bam_record.sequence().iter().collect();
                if seq.is_empty() {
                    warn(format!(
                        "skipping read '{name}' in '{}': it has no sequence",
                        self.filename.display()
                    ));
                    self.skipped += 1;
                    continue;
                }
                let seq = seq.to_ascii_uppercase();
                // A reverse-complemented flag on an unaligned record means the stored
                // sequence was flipped; put it back as it was sequenced.
                let seq = if flags.is_reverse_complemented() {
                    dna::revcomp(&seq)
                } else {
                    seq
                };
                return Ok(Some(RawRead { name, seq }));
            },
        }
    }
}

/// Summarise records dropped while loading, so warnings scrolling past are not the only
/// record of how many reads actually made it through.
fn report_skipped(skipped: usize, loaded: usize, filename: &Path) {
    if skipped > 0 {
        warn(format!(
            "skipped {skipped} unusable record(s) in '{}'; {loaded} read(s) loaded",
            filename.display()
        ));
    }
}

/// Load every read from a sequences file at once.
///
/// `run` streams the file in chunks instead; this is for callers that want the whole
/// set in memory, and is what the tests exercise.
#[cfg(test)]
fn load_reads(filename: &Path) -> Result<LoadedReads> {
    let mut stream = ReadStream::open(filename)?;
    let mut reads = Vec::new();
    loop {
        let chunk = stream.next_chunk(usize::MAX, usize::MAX)?;
        if chunk.is_empty() {
            break;
        }
        reads.extend(chunk);
    }
    report_skipped(stream.skipped, reads.len(), filename);
    Ok(LoadedReads {
        skipped: stream.skipped,
        reads,
    })
}

/// Cut out the relevant region of a read and if necessary generate the reverse complement
/// so that a sequence from 'start' to 'end' segments is always returned.
fn cut_reorient_seq(seq: &[u8], start_idx: usize, end_idx: usize, revcomp: bool) -> Vec<u8> {
    let new_seq: Vec<u8>;
    if start_idx < end_idx {
        new_seq = seq[start_idx..end_idx].to_vec();
    } else {
        if revcomp {
            let rc_seq = dna::revcomp(seq);
            new_seq = rc_seq[(seq.len() - start_idx)..(seq.len() - end_idx)].to_vec();
        } else {
            let mut rev_seq = seq.to_vec();
            rev_seq.reverse();
            new_seq = rev_seq[(seq.len() - start_idx)..(seq.len() - end_idx)].to_vec();
        }
    }
    new_seq
}

/// Scoring for the start/end anchors: +2 match, -1 mismatch, and an affine gap where a
/// run of k gaps costs -2 - 2k.
const ANCHOR_SCORING: (i32, i32, i32, i32) = (2, -1, -2, -2);

/// A semi-global alignment reduced to what the caller actually needs: its score and the
/// span of the read it covers.
pub struct AnchorAlignment {
    pub score: i32,
    pub ystart: usize,
    pub yend: usize,
}

/// Best semi-global alignment of `query` within `target`, in space linear in the target.
///
/// Mirrors the affine gap model of the library aligner it replaces, where a run of `k`
/// gaps costs `gap_open + gap_extend * k` - so with both set to -2 a single gap costs -4,
/// not -2. Three running rows track alignments ending in a match, in a gap in the target,
/// and in a gap in the query, each carrying the target position its path began at.
fn best_semiglobal(query: &[u8], target: &[u8], scoring: (i32, i32, i32, i32)) -> AnchorAlignment {
    let (match_score, mismatch_score, gap_open, gap_extend) = scoring;
    let query_len = query.len();
    let target_len = target.len();
    const NEG: i32 = i32::MIN / 4;

    // m: ends aligned to a base. i: ends in a gap in the target (query consumed).
    // d: ends in a gap in the query (target consumed). `*_s` are the matching starts.
    let mut m = vec![0i32; target_len + 1];
    let mut i_row = vec![NEG; target_len + 1];
    let mut d = vec![NEG; target_len + 1];
    let mut m_s: Vec<usize> = (0..=target_len).collect();
    let mut i_s = vec![0usize; target_len + 1];
    let mut d_s = vec![0usize; target_len + 1];
    let (mut nm, mut ni, mut nd) = (
        vec![0i32; target_len + 1],
        vec![NEG; target_len + 1],
        vec![NEG; target_len + 1],
    );
    let (mut nm_s, mut ni_s, mut nd_s) = (
        vec![0usize; target_len + 1],
        vec![0usize; target_len + 1],
        vec![0usize; target_len + 1],
    );

    for qi in 1..=query_len {
        nm[0] = NEG;
        ni[0] = gap_open + gap_extend * qi as i32;
        nd[0] = NEG;
        nm_s[0] = 0;
        ni_s[0] = 0;
        nd_s[0] = 0;
        let query_char = query[qi - 1];
        for j in 1..=target_len {
            let sub = if query_char == target[j - 1] {
                match_score
            } else {
                mismatch_score
            };
            // Diagonal: whichever predecessor state scored best, preferring a match.
            let (best_prev, best_start) = pick3(
                (m[j - 1], m_s[j - 1]),
                (i_row[j - 1], i_s[j - 1]),
                (d[j - 1], d_s[j - 1]),
            );
            nm[j] = best_prev.saturating_add(sub);
            nm_s[j] = best_start;
            // Gap in the target: consume a query base, staying in the same column.
            let open_i = m[j].saturating_add(gap_open + gap_extend);
            let extend_i = i_row[j].saturating_add(gap_extend);
            if open_i >= extend_i {
                ni[j] = open_i;
                ni_s[j] = m_s[j];
            } else {
                ni[j] = extend_i;
                ni_s[j] = i_s[j];
            }
            // Gap in the query: consume a target base, staying in the same row.
            let open_d = nm[j - 1].saturating_add(gap_open + gap_extend);
            let extend_d = nd[j - 1].saturating_add(gap_extend);
            if open_d >= extend_d {
                nd[j] = open_d;
                nd_s[j] = nm_s[j - 1];
            } else {
                nd[j] = extend_d;
                nd_s[j] = nd_s[j - 1];
            }
        }
        std::mem::swap(&mut m, &mut nm);
        std::mem::swap(&mut i_row, &mut ni);
        std::mem::swap(&mut d, &mut nd);
        std::mem::swap(&mut m_s, &mut nm_s);
        std::mem::swap(&mut i_s, &mut ni_s);
        std::mem::swap(&mut d_s, &mut nd_s);
    }

    // The whole query must be consumed; trailing target is free, so scan every end.
    let mut best = AnchorAlignment {
        score: NEG,
        ystart: 0,
        yend: 0,
    };
    for j in 0..=target_len {
        let (score, start) = if m[j] >= i_row[j] {
            (m[j], m_s[j])
        } else {
            (i_row[j], i_s[j])
        };
        if score > best.score {
            best = AnchorAlignment {
                score,
                ystart: start,
                yend: j,
            };
        }
    }
    best
}

/// Best of three (score, start) candidates, preferring the earlier on a tie.
fn pick3(a: (i32, usize), b: (i32, usize), c: (i32, usize)) -> (i32, usize) {
    let ab = if a.0 >= b.0 { a } else { b };
    if ab.0 >= c.0 { ab } else { c }
}

/// The library aligner this replaced, kept only so the tests can check the two still
/// agree. Its affine gap model is what `best_semiglobal` reproduces.
#[cfg(test)]
fn align_segment(segment_seq: &[u8], read_seq: &[u8]) -> bio::alignment::Alignment {
    use bio::alignment::pairwise::Aligner;
    let score = |a: u8, b: u8| if a == b { 2i32 } else { -1i32 };
    // gap open score: -2, gap extension score: -2
    let mut aligner = Aligner::with_capacity(segment_seq.len(), read_seq.len(), -2, -2, &score);
    aligner.semiglobal(segment_seq, read_seq)
}

/// The 'start' and 'end' anchors with their reverse complements, computed once per run.
type AnchorSeqs = (Vec<u8>, Vec<u8>, Vec<u8>, Vec<u8>);

/// Work out why an anchored read was rejected.
///
/// Each anchor is judged by its better score across both strands rather than by the
/// branch taken, since the branch is chosen before the score check and would blame both
/// anchors whenever either was missing. Both present but still failing means they were
/// found facing opposite ways.
fn anchor_failure(start_found: bool, end_found: bool) -> Rejected {
    match (start_found, end_found) {
        (false, false) => Rejected::NeitherAnchorFound,
        (false, true) => Rejected::StartAnchorNotFound,
        (true, false) => Rejected::EndAnchorNotFound,
        (true, true) => Rejected::AnchorsOnOppositeStrands,
    }
}

/// Classify which segments and their orientations are present within a read.
///
/// With `start_end_segs` the read has already been trimmed to run from 'start' to 'end',
/// so those two bound the segment string instead of being searched for within it.
/// Without it the read is classified as given and a segment named 'start' or 'end' is
/// treated like any other.
fn classify_read_segments(
    segments: &HashMap<String, Vec<u8>>,
    read_seq: &[u8],
    min_norm_score: f32,
    start_end_segs: bool,
) -> String {
    let aligner = MultiTracebackAligner::new(
        2,  // match score
        -1, // mismatch score
        -2, // gap score
    );
    let mut all_alignments: Vec<MultiAlignment> = Vec::new();
    for (seg_name, seg_seq) in segments {
        if start_end_segs && (seg_name == "start" || seg_name == "end") {
            continue;
        }
        // Rejected when the segments file is loaded, but guard anyway: a zero-length
        // segment scores zero everywhere, so it would "match" at every position.
        if seg_seq.is_empty() {
            continue;
        }
        let min_score = (seg_seq.len() as f32 * min_norm_score) as i32;
        all_alignments.append(&mut aligner.align_multiple(seg_name, seg_seq, read_seq, min_score));
        all_alignments.append(&mut aligner.align_multiple(
            &format!("{}*", seg_name),
            &dna::revcomp(seg_seq),
            read_seq,
            min_score,
        ));
    }
    // Highest normalised score first; total_cmp orders every float, so a stray NaN can
    // never panic here. Position and name break ties: segments live in a HashMap whose
    // iteration order varies between runs and the sort is stable, so without an explicit
    // tiebreak the greedy filtering below could drop a different one each run.
    all_alignments.sort_by(|a, b| {
        b.norm_score
            .total_cmp(&a.norm_score)
            .then_with(|| a.target_start.cmp(&b.target_start))
            .then_with(|| a.target_end.cmp(&b.target_end))
            .then_with(|| a.name.cmp(&b.name))
    });
    // Accept greedily, best first. Segments in real constructs often abut or share a
    // few bases, so a small overlap is tolerated; beyond that the weaker hit is dropped.
    let allowed_overlap: usize = 3;
    if all_alignments.len() > 1 {
        for idx1 in 0..(all_alignments.len() - 1) {
            if all_alignments[idx1].filter {
                continue;
            } else {
                for idx2 in (idx1 + 1)..all_alignments.len() {
                    if all_alignments[idx2].filter {
                        continue;
                    } else {
                        let s1 = all_alignments[idx1].target_start;
                        let s2 = all_alignments[idx2].target_start;
                        let e1 = all_alignments[idx1].target_end;
                        let e2 = all_alignments[idx2].target_end;
                        if (s2 < s1 && e2 <= s1) || (s2 >= e1 && e2 > e1) {
                            continue;
                        }
                        if s2 <= s1 && e2 >= e1 {
                            all_alignments[idx2].filter = true;
                            continue;
                        }
                        if s2 >= s1 && s2 <= e1 && e1 - s2 > allowed_overlap {
                            all_alignments[idx2].filter = true;
                            continue;
                        }
                        if e2 >= s1 && e2 <= e1 && e2 - s1 > allowed_overlap {
                            all_alignments[idx2].filter = true;
                            continue;
                        }
                    }
                }
            }
        }
    }
    all_alignments.sort_by_key(|a| a.target_start);
    // Generate the segment string, bounded by the anchors only if the user supplied them
    let mut seg_names: Vec<&str> = Vec::new();
    if start_end_segs {
        seg_names.push("start");
    }
    for aln in &all_alignments {
        if aln.filter {
            continue;
        }
        seg_names.push(aln.name.as_str());
    }
    if start_end_segs {
        seg_names.push("end");
    }
    seg_names.join("-")
}

/// Classify segments across all reads loaded from a sequences file (FASTQ or BAM).
fn process_reads(
    records: Vec<RawRead>,
    segments: &HashMap<String, Vec<u8>>,
    len_check: (usize, usize),
    min_norm_score: f32,
    circular: bool,
    start_end_segs: bool,
) -> Result<(Vec<SegmentedRead>, RunSummary)> {
    // Precompute start/end revcomps once rather than per-read
    let start_end_seqs: Option<AnchorSeqs> = if start_end_segs {
        let anchor = |name: &str| -> Result<Vec<u8>> {
            segments.get(name).cloned().ok_or_else(|| {
                format!(
                    "--start-end-segs was given, but the segments file has no segment named \
                     '{name}'. Add it, or drop --start-end-segs to classify reads as they are."
                )
            })
        };
        let start_seq = anchor("start")?;
        let end_seq = anchor("end")?;
        let start_rc_seq = dna::revcomp(&start_seq);
        let end_rc_seq = dna::revcomp(&end_seq);
        Some((start_seq, start_rc_seq, end_seq, end_rc_seq))
    } else {
        None
    };

    let clean_seqs: Vec<Outcome> = records
        .par_iter()
        .map(|record| {
            let read_name = record.name.clone();
            let read_seq = record.seq.clone();

            let read_len = read_seq.len();
            if len_check != (0, 0) && read_len < len_check.0 {
                return Outcome::Rejected(Rejected::TooShort);
            }
            if len_check != (0, 0) && read_len > len_check.1 {
                return Outcome::Rejected(Rejected::TooLong);
            }

            if !start_end_segs {
                // No anchors: classify the read as it arrived, on whichever strand.
                let clean_seq = cut_reorient_seq(&read_seq, 0, read_seq.len(), false);
                let seg_str =
                    classify_read_segments(segments, &clean_seq, min_norm_score, start_end_segs);
                Outcome::Classified(SegmentedRead {
                    name: read_name,
                    segments: seg_str,
                })
            } else if let Some((start_seq, start_rc_seq, end_seq, end_rc_seq)) =
                start_end_seqs.as_ref()
            {
                let anchor = |seq: &[u8]| best_semiglobal(seq, &read_seq, ANCHOR_SCORING);
                let align_start = anchor(start_seq);
                let align_end = anchor(end_seq);
                let align_start_rc = anchor(start_rc_seq);
                let align_end_rc = anchor(end_rc_seq);

                let start_min = (start_seq.len() as f32 * min_norm_score) as i32;
                let end_min = (end_seq.len() as f32 * min_norm_score) as i32;
                if align_start.score > align_start_rc.score && align_end.score > align_end_rc.score
                {
                    // Read in the correct orientation
                    if align_start.score >= start_min && align_end.score >= end_min {
                        let new_read_seq: Vec<u8>;
                        let start_seg_pos: usize;
                        let end_seg_pos: usize;
                        if circular && align_end.ystart < align_start.ystart {
                            let p1 = &read_seq[align_start.ystart..];
                            let p2 = &read_seq[0..align_end.yend];
                            new_read_seq = [p1, p2].concat();
                            start_seg_pos = 0;
                            end_seg_pos = new_read_seq.len();
                        } else if align_start.ystart < align_end.yend {
                            new_read_seq = read_seq.clone();
                            start_seg_pos = align_start.ystart;
                            end_seg_pos = align_end.yend;
                        } else {
                            // 'start' does not precede 'end' and this isn't a circular
                            // wrap-around: the read is inconsistent, reject it.
                            return Outcome::Rejected(Rejected::AnchorsOutOfOrder);
                        }
                        let clean_seq =
                            cut_reorient_seq(&new_read_seq, start_seg_pos, end_seg_pos, false);
                        let seg_str = classify_read_segments(
                            segments,
                            &clean_seq,
                            min_norm_score,
                            start_end_segs,
                        );
                        Outcome::Classified(SegmentedRead {
                            name: read_name,
                            segments: seg_str,
                        })
                    } else {
                        Outcome::Rejected(anchor_failure(
                            align_start.score.max(align_start_rc.score) >= start_min,
                            align_end.score.max(align_end_rc.score) >= end_min,
                        ))
                    }
                } else {
                    // Read in the reverse orientation
                    if align_start_rc.score >= start_min && align_end_rc.score >= end_min {
                        let new_read_seq: Vec<u8>;
                        let start_seg_pos: usize;
                        let end_seg_pos: usize;
                        if circular && align_start_rc.ystart < align_end_rc.ystart {
                            let p1 = &read_seq[align_end_rc.ystart..];
                            let p2 = &read_seq[0..align_start_rc.yend];
                            new_read_seq = [p1, p2].concat();
                            end_seg_pos = 0;
                            start_seg_pos = new_read_seq.len();
                        } else if align_start_rc.yend > align_end_rc.ystart {
                            new_read_seq = read_seq.clone();
                            start_seg_pos = align_start_rc.yend;
                            end_seg_pos = align_end_rc.ystart;
                        } else {
                            // 'end' does not precede 'start' and this isn't a circular
                            // wrap-around: the read is inconsistent, reject it.
                            return Outcome::Rejected(Rejected::AnchorsOutOfOrder);
                        }
                        let clean_seq =
                            cut_reorient_seq(&new_read_seq, start_seg_pos, end_seg_pos, true);
                        let seg_str = classify_read_segments(
                            segments,
                            &clean_seq,
                            min_norm_score,
                            start_end_segs,
                        );
                        Outcome::Classified(SegmentedRead {
                            name: read_name,
                            segments: seg_str,
                        })
                    } else {
                        Outcome::Rejected(anchor_failure(
                            align_start.score.max(align_start_rc.score) >= start_min,
                            align_end.score.max(align_end_rc.score) >= end_min,
                        ))
                    }
                }
            } else {
                Outcome::Rejected(Rejected::NeitherAnchorFound)
            }
        })
        .collect();

    let mut summary = RunSummary::default();
    let mut classified = Vec::new();
    for outcome in clean_seqs {
        match outcome {
            Outcome::Classified(read) => {
                summary.classified += 1;
                classified.push(read);
            }
            Outcome::Rejected(reason) => summary.count(reason),
        }
    }
    Ok((classified, summary))
}

/// One alignment of a segment within a read: where it sits, how well it scored, and
/// whether overlap resolution has discarded it.
#[derive(Debug, Clone)]
pub struct MultiAlignment {
    pub name: String,
    pub score: i32,
    pub norm_score: f32,
    pub target_start: usize,
    pub target_end: usize,
    pub filter: bool,
}

/// Semi-global aligner that reports every place a segment matches a read, not just
/// the best one.
pub struct MultiTracebackAligner {
    match_score: i32,
    mismatch_score: i32,
    gap_score: i32,
}

impl MultiTracebackAligner {
    /// Build an aligner scoring matches, mismatches and gaps as given.
    pub fn new(match_score: i32, mismatch_score: i32, gap_score: i32) -> Self {
        Self {
            match_score,
            mismatch_score,
            gap_score,
        }
    }

    /// Find every semi-global alignment of `query` within `target` scoring at least
    /// `min_score`, one per end position.
    ///
    /// Runs in space linear in the target: rather than keep the whole (query x target)
    /// matrix only to trace back where each alignment began, every cell carries that
    /// start position forward, so it falls out of the recurrence directly.
    ///
    /// Ties break diagonal, then up, then left - the precedence the old backward
    /// traceback used, and what makes this exactly equivalent to it rather than merely
    /// similar. `src/tests.rs` keeps that original as a reference oracle.
    pub fn align_multiple(
        &self,
        seg_name: &str,
        query: &[u8],
        target: &[u8],
        min_score: i32,
    ) -> Vec<MultiAlignment> {
        let query_len = query.len();
        let target_len = target.len();

        // Row 0. Leading gaps in the target are free, so every column is a possible
        // start and scores zero.
        let mut prev_score: Vec<i32> = vec![0; target_len + 1];
        let mut prev_start: Vec<usize> = (0..=target_len).collect();
        let mut curr_score: Vec<i32> = vec![0; target_len + 1];
        let mut curr_start: Vec<usize> = vec![0; target_len + 1];

        for i in 1..=query_len {
            // Column 0: the query consumed so far is gapped against nothing.
            curr_score[0] = self.gap_score * i as i32;
            curr_start[0] = 0;
            let query_char = query[i - 1];
            for j in 1..=target_len {
                let match_score = if query_char == target[j - 1] {
                    self.match_score
                } else {
                    self.mismatch_score
                };
                let diagonal_score = prev_score[j - 1] + match_score;
                let up_score = prev_score[j] + self.gap_score;
                let left_score = curr_score[j - 1] + self.gap_score;
                let max_score = diagonal_score.max(up_score).max(left_score);
                curr_score[j] = max_score;
                curr_start[j] = if diagonal_score == max_score {
                    prev_start[j - 1]
                } else if up_score == max_score {
                    prev_start[j]
                } else {
                    curr_start[j - 1]
                };
            }
            std::mem::swap(&mut prev_score, &mut curr_score);
            std::mem::swap(&mut prev_start, &mut curr_start);
        }

        // After the loop `prev_*` holds the last row computed, i.e. row `query_len`.
        // With an empty query that is row 0, which is what the caller should see.
        let mut alignments: Vec<MultiAlignment> = Vec::new();
        for j in 0..=target_len {
            let score = prev_score[j];
            if score >= min_score {
                alignments.push(MultiAlignment {
                    name: seg_name.to_string(),
                    score,
                    norm_score: score as f32 / query_len as f32,
                    target_start: prev_start[j],
                    target_end: j,
                    filter: false,
                });
            }
        }
        alignments
    }
}

#[derive(Parser)]
#[command(name = "segment")]
#[command(about = "Segment alignment tool")]
#[command(version)]
/// Command line arguments.
struct Cli {
    /// File with segment sequences to search for (FASTA format, one record per segment)
    #[arg(long)]
    segments: String,

    /// File with reads to process: FASTQ, gzipped FASTQ, or unaligned BAM. The format is
    /// detected from the file contents, so the extension does not matter
    #[arg(long)]
    sequences: String,

    /// Minimum read length to process
    #[arg(long, default_value_t = 100)]
    min_seq_len: usize,

    /// Maximum read length to process
    #[arg(long, default_value_t = 100000)]
    max_seq_len: usize,

    /// Use the 'start' and 'end' segments to anchor, trim and orient each read. They
    /// then bound the segment string rather than being reported within it. Without
    /// this, reads are classified as given and any segment named 'start' or 'end' is
    /// treated as an ordinary segment
    #[arg(short = 'd', long, default_value_t = false)]
    start_end_segs: bool,

    /// Minimum normalised score, between 0.0 and 2.0 (perfect match)
    #[arg(short = 's', long, default_value_t = 1.5)]
    min_norm_score: f32,

    /// Sequences to classify are circular (i.e., from plasmids)
    #[arg(short = 'c', long, default_value_t = false)]
    circular: bool,

    /// Threads to use
    #[arg(short = 't', long, default_value_t = 1)]
    threads: usize,

    /// Output file
    #[arg(short = 'o', long, default_value = "output.txt")]
    output: String,
}

/// Run the tool, reporting any failure on stderr and exiting non-zero.
fn main() {
    if let Err(message) = run() {
        eprintln!("error: {message}");
        std::process::exit(1);
    }
}

/// Load the inputs, classify every read, write the results and report the summary.
fn run() -> Result<()> {
    let cli = Cli::parse();

    if !(0.0..=2.0).contains(&cli.min_norm_score) {
        return Err(format!(
            "--min-norm-score must be between 0.0 and 2.0 (2.0 is a perfect match), but {} was given.",
            cli.min_norm_score
        ));
    }
    if cli.min_seq_len > cli.max_seq_len {
        return Err(format!(
            "--min-seq-len ({}) is greater than --max-seq-len ({}), so no read could ever be kept.",
            cli.min_seq_len, cli.max_seq_len
        ));
    }

    rayon::ThreadPoolBuilder::new()
        .num_threads(cli.threads)
        .build_global()
        .map_err(|e| {
            format!(
                "Could not start a thread pool of {} thread(s): {e}",
                cli.threads
            )
        })?;

    let segments = load_fasta(Path::new(&cli.segments))?;
    let sequences = Path::new(&cli.sequences);
    let mut stream = ReadStream::open(sequences)?;

    let f = File::create(Path::new(&cli.output))
        .map_err(|e| format!("Could not create output file '{}': {e}", cli.output))?;
    let mut f = BufWriter::new(f);

    // Progress is measured in bytes of the file consumed. That works the same whether the
    // input is plain, gzipped or BAM, and needs no counting pass over the file first.
    let total_bytes = std::fs::metadata(sequences).map(|m| m.len()).unwrap_or(0);
    let progress = ProgressBar::new(total_bytes);
    progress.set_style(
        ProgressStyle::with_template(
            "{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {bytes}/{total_bytes} ({percent}%, {bytes_per_sec}, ETA {eta}) {msg}",
        )
        .unwrap()
        .progress_chars("#>-"),
    );

    // Read, classify and write one bounded chunk at a time. Chunks are handled in order
    // and rayon preserves order within a chunk, so the output is byte for byte what
    // processing the whole file at once would produce.
    let mut summary = RunSummary::default();
    let mut reads_seen = 0usize;
    loop {
        let chunk = stream.next_chunk(MAX_CHUNK_BASES, MAX_CHUNK_READS)?;
        if chunk.is_empty() {
            break;
        }
        reads_seen += chunk.len();
        let (clean_seqs, chunk_summary) = process_reads(
            chunk,
            &segments,
            (cli.min_seq_len, cli.max_seq_len),
            cli.min_norm_score,
            cli.circular,
            cli.start_end_segs,
        )?;
        for r in clean_seqs {
            writeln!(f, "{}\t{}", r.name, r.segments)
                .map_err(|e| format!("Could not write to output file '{}': {e}", cli.output))?;
        }
        summary.merge(chunk_summary);
        // Advanced once the chunk is written, so the bar means "done", not "read ahead".
        progress.set_position(stream.bytes_read().min(total_bytes));
        progress.set_message(format!("{reads_seen} reads"));
    }
    progress.set_position(total_bytes);
    progress.finish_and_clear();

    f.flush()
        .map_err(|e| format!("Could not finish writing output file '{}': {e}", cli.output))?;

    report_skipped(stream.skipped, reads_seen, sequences);
    if reads_seen == 0 {
        warn(format!(
            "no usable reads were loaded from '{}'",
            cli.sequences
        ));
    }
    summary.unreadable = stream.skipped;
    eprint!("{}", summary.render());
    Ok(())
}
