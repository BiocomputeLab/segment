// Segment alignment tool
// Author: Thomas E. Gorochowski <tom@chofski.co.uk>

use bio::alphabets::dna;
use bio::io::fasta::FastaRead;
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
use std::time::{Duration, Instant, SystemTime};

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

/// Nucleotide sets, one bit per base. An IUPAC ambiguity code is the union of the bases
/// it stands for, so two bases match exactly when their sets intersect. On plain ACGT
/// this is byte equality, which is what the aligners did before ambiguity codes were
/// understood at all.
const BASE_A: u8 = 0b0001;
const BASE_C: u8 = 0b0010;
const BASE_G: u8 = 0b0100;
const BASE_T: u8 = 0b1000;
/// Every base: what `N` permits, and the widest set any code can name.
const BASE_ANY: u8 = BASE_A | BASE_C | BASE_G | BASE_T;

/// Every IUPAC nucleotide code and the bases it permits. `U` is folded into `T` when
/// segments are loaded, so it never reaches this table.
const IUPAC_CODES: [(u8, u8); 15] = [
    (b'A', BASE_A),
    (b'C', BASE_C),
    (b'G', BASE_G),
    (b'T', BASE_T),
    (b'R', BASE_A | BASE_G),
    (b'Y', BASE_C | BASE_T),
    (b'S', BASE_C | BASE_G),
    (b'W', BASE_A | BASE_T),
    (b'K', BASE_G | BASE_T),
    (b'M', BASE_A | BASE_C),
    (b'B', BASE_C | BASE_G | BASE_T),
    (b'D', BASE_A | BASE_G | BASE_T),
    (b'H', BASE_A | BASE_C | BASE_T),
    (b'V', BASE_A | BASE_C | BASE_G),
    (b'N', BASE_ANY),
];

/// Byte to nucleotide set. Anything that is not an IUPAC code maps to the empty set and
/// so matches nothing - not even a copy of itself. Segments are rejected on load if they
/// contain one; reads are not, because a stray character in noisy data should cost that
/// position and nothing more.
static BASE_MASK: [u8; 256] = build_base_masks();

/// Build [`BASE_MASK`] at compile time, accepting either case.
const fn build_base_masks() -> [u8; 256] {
    let mut masks = [0u8; 256];
    let mut i = 0;
    while i < IUPAC_CODES.len() {
        let code = IUPAC_CODES[i].0;
        let mask = IUPAC_CODES[i].1;
        masks[code as usize] = mask;
        masks[(code + 32) as usize] = mask; // the lower-case spelling
        i += 1;
    }
    masks
}

/// Encode a sequence as one nucleotide set per base, ready for the aligners.
fn base_masks(seq: &[u8]) -> Vec<u8> {
    seq.iter().map(|&base| BASE_MASK[base as usize]).collect()
}

/// A segment to search for: the sequence, and the normalised score a hit has to reach to
/// count. The threshold is the global `--min-norm-score` unless the segments file named
/// a different one for this segment, so nothing downstream has to know which it was.
#[derive(Debug, Clone)]
struct Segment {
    seq: Vec<u8>,
    min_norm_score: f32,
}

/// How to write a per-segment score, quoted in most of the messages below.
const SCORE_SYNTAX: &str = "With --per-segment-scores the square brackets hold that segment's \
     minimum normalised score, between 0.0 and 2.0, as in 'SEG1[1.5]'.";

/// Split a FASTA record ID into the segment name and the threshold in square brackets
/// after it: `SEG1[1.5]` is the segment `SEG1`, found at 1.5.
///
/// A trailing bracketed value is always taken off the name, so a segment is reported
/// under the same name however the run is configured. `--per-segment-scores` decides
/// only whether the value inside is *used*: with the flag it must be a valid score and
/// becomes that segment's threshold, without it the value is ignored and the segment is
/// held to the global `--min-norm-score` like any other.
///
/// With the flag on, a bracket that cannot be read as a score is a mistake worth stopping
/// for rather than something to quietly fold into the name: a segment silently renamed is
/// a segment silently never found, and nothing later in the run would hint at why it went
/// missing.
///
/// Errors are the detail only. `load_fasta` prefixes them with the record they came from.
fn parse_segment_name(id: &str, per_segment_scores: bool) -> Result<(String, Option<f32>)> {
    let Some((name, written)) = id.strip_suffix(']').and_then(|rest| rest.rsplit_once('[')) else {
        // No closed bracket group at the end. With the flag on, a bracket left open is a
        // typo rather than a name: saying so beats accepting 'SEG1[1.5' as a segment
        // nothing will ever be found under. Without it there is no syntax to get wrong.
        if per_segment_scores && let Some(open) = id.rfind('[') {
            return Err(if id[open..].contains(']') {
                format!(
                    "the score has to come last, but '{}' follows the closing bracket. {}",
                    &id[id.rfind(']').unwrap_or(open) + 1..],
                    SCORE_SYNTAX
                )
            } else {
                format!("the square bracket is never closed. {SCORE_SYNTAX}")
            });
        }
        return Ok((id.to_string(), None));
    };
    // Broken whether or not the brackets are being read: stripping them would leave a
    // segment with no name at all.
    if name.is_empty() {
        return Err(
            "there is no name in front of the brackets. Write the name first, as in 'SEG1[1.5]'."
                .to_string(),
        );
    }
    if !per_segment_scores {
        // The brackets come off regardless, but nothing inside them is interpreted, so
        // whatever they hold is neither validated nor applied.
        return Ok((name.to_string(), None));
    }
    // Reported back as written rather than as the reparsed float, so a message quotes
    // what is actually in the file: '2.50' does not come back as '2.5'.
    let written = written.trim();
    if written.is_empty() {
        return Err(format!("the square brackets are empty. {SCORE_SYNTAX}"));
    }
    let score: f32 = written
        .parse()
        .map_err(|_| format!("'{written}' is not a number. {SCORE_SYNTAX}"))?;
    // Caught before the range check below, whose advice would otherwise be nonsense: a
    // NaN is neither too high nor too low, it is simply not a score.
    if !score.is_finite() {
        return Err(format!("'{written}' is not a usable score. {SCORE_SYNTAX}"));
    }
    if !(0.0..=2.0).contains(&score) {
        // Which end it fell off says what would have happened had it been allowed, which
        // is more use than repeating the bounds back.
        let consequence = if score > 2.0 {
            "no alignment can score above 2.0, so the segment could never be found"
        } else {
            "every alignment already scores at least 0.0, so the segment would be found \
             everywhere"
        };
        return Err(format!(
            "a minimum normalised score of '{written}' is outside the range 0.0 to 2.0, where 2.0 \
             is a perfect match - {consequence}."
        ));
    }
    Ok((name.to_string(), Some(score)))
}

/// Whether a record ID looks like it was meant to name a score, for warning about a
/// segments file that uses the syntax without the flag that turns it on.
fn looks_like_a_scored_name(id: &str) -> bool {
    id.strip_suffix(']')
        .and_then(|rest| rest.rsplit_once('['))
        .is_some_and(|(name, score)| !name.is_empty() && score.trim().parse::<f32>().is_ok())
}

/// Loads a FASTA file containing the `start` and `end` sequences plus and
/// other segment sequences that should be used when classifying a read.
///
/// Sequences are upper-cased on the way in so that a soft-masked or lower-case FASTA
/// still matches reads, and `U` is folded into `T` so an RNA spelling works too. Every
/// segment leaves here with a threshold already resolved - its own if the file named
/// one, `default_min_norm_score` otherwise - so no later stage needs the global value or
/// has to know which segments overrode it.
///
/// Every problem here is fatal: a segment that cannot be used is a mistake in the input,
/// and carrying on would silently change what gets reported.
fn load_fasta(
    filename: &Path,
    default_min_norm_score: f32,
    per_segment_scores: bool,
) -> Result<HashMap<String, Segment>> {
    let mut fasta_data: HashMap<String, Segment> = HashMap::new();
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
        // The name is checked for uniqueness after the score is stripped off it, so
        // 'SEG1[1.5]' and 'SEG1[1.7]' are the same segment asking for two thresholds.
        //
        // Failures are reported against the raw ID and the record number rather than the
        // parsed name, since with a malformed bracket there may not be a usable name yet,
        // and the record number is what makes it findable in a file of hundreds.
        let (seg_name, seg_score) =
            parse_segment_name(record.id(), per_segment_scores).map_err(|detail| {
                format!(
                    "Segment '{}' in '{}' (record {}): {detail}",
                    record.id(),
                    filename.display(),
                    index + 1,
                )
            })?;
        // Someone who writes the syntax but forgets the flag gets the threshold they
        // asked for silently ignored, with nothing in the run to say so. Cheap to spot.
        if !per_segment_scores && looks_like_a_scored_name(record.id()) {
            warn(format!(
                "segment '{}' in '{}' names a minimum normalised score, but --per-segment-scores \
                 was not given, so it is ignored and the segment is held to --min-norm-score \
                 like any other. The brackets are dropped from its name either way.",
                record.id(),
                filename.display(),
            ));
        }
        let mut seg_seq = record.seq().to_ascii_uppercase();
        // Done here rather than in the mask table so that reverse complementing stays
        // correct: `U` has no complement of its own, and would be left as it is.
        for base in seg_seq.iter_mut() {
            if *base == b'U' {
                *base = b'T';
            }
        }
        if seg_seq.is_empty() {
            return Err(format!(
                "Segment '{seg_name}' in '{}' has no sequence. Every segment needs at least one base.",
                filename.display()
            ));
        }
        // A character outside the alphabet matches nothing, so leaving it in place would
        // quietly stop the segment ever being found. Say so instead.
        if let Some(offset) = seg_seq.iter().position(|&b| BASE_MASK[b as usize] == 0) {
            return Err(format!(
                "Segment '{seg_name}' in '{}' has '{}' at position {}, which is not a nucleotide. \
                 Segments may use A, C, G and T, the IUPAC ambiguity codes R, Y, S, W, K, M, B, D, \
                 H, V and N, or U in place of T.",
                filename.display(),
                seg_seq[offset].escape_ascii(),
                offset + 1,
            ));
        }
        // Every base ambiguous means the segment matches perfectly at every position of
        // every read, which is never what anyone meant to ask for.
        if seg_seq.iter().all(|&b| BASE_MASK[b as usize] == BASE_ANY) {
            return Err(format!(
                "Segment '{seg_name}' in '{}' is entirely ambiguous, so it would match at every \
                 position of every read. Give it at least one base that is not N.",
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
        fasta_data.insert(
            seg_name,
            Segment {
                seq: seg_seq,
                min_norm_score: seg_score.unwrap_or(default_min_norm_score),
            },
        );
    }
    if fasta_data.is_empty() {
        return Err(format!(
            "No segments found in '{}'. Expected a FASTA file with one record per segment.",
            filename.display()
        ));
    }
    Ok(fasta_data)
}

/// Segments ambiguous enough to be found in sequence that does not contain them, each
/// paired with the fraction of random bases it matches and the fraction its own
/// threshold demands. Sorted by name.
///
/// An ungapped alignment reaching a normalised score of `s` needs a fraction
/// `(s + 1) / 3` of its bases to match, since a match scores +2 and a mismatch -1. An
/// ambiguity code matches a random base with probability `|set| / 4`, so when a
/// segment's average over its bases reaches what the threshold demands, the segment
/// clears that threshold on sequence picked at random. That is not an error - a mostly
/// degenerate probe is a legitimate thing to search for - but the results will be full
/// of it, so it should not come as a surprise.
fn degenerate_segments(segments: &HashMap<String, Segment>) -> Vec<(String, f32, f32)> {
    // Segments live in a HashMap whose iteration order varies between runs; sort so the
    // warnings always come out in the same order.
    let mut names: Vec<&String> = segments.keys().collect();
    names.sort();
    let mut degenerate = Vec::new();
    for name in names {
        let segment = &segments[name];
        // Rejected when the segments file is loaded, but guard anyway: dividing by a
        // zero length below would give a NaN that compares false against everything.
        if segment.seq.is_empty() {
            continue;
        }
        // Judged against this segment's own threshold, so a segment that raised its own
        // to cover its ambiguity is not reported.
        let required = required_identity(segment.min_norm_score);
        let chance = chance_match_fraction(&segment.seq);
        if chance >= required {
            degenerate.push((name.clone(), chance, required));
        }
    }
    degenerate
}

/// The fraction of an ungapped segment's bases that must match to reach `min_norm_score`.
fn required_identity(min_norm_score: f32) -> f32 {
    (min_norm_score + 1.0) / 3.0
}

/// How much of a random sequence a segment matches, averaged over its bases. A spelled
/// out base matches one base in four; an ambiguity code matches as many as it permits.
fn chance_match_fraction(seq: &[u8]) -> f32 {
    seq.iter()
        .map(|&base| BASE_MASK[base as usize].count_ones() as f32 / 4.0)
        .sum::<f32>()
        / seq.len() as f32
}

/// Structure to hold the key sequence and segment information about a read.
struct SegmentedRead {
    name: String,
    segments: String,
    /// Filled only with `--detailed-output`, which is what decides whether the extra
    /// columns are written. The extracted sequence roughly doubles what a chunk holds
    /// while it is being written, so it is not carried when nothing will print it.
    detail: Option<ReadDetail>,
}

/// The extra columns `--detailed-output` adds: where each segment sits in the extracted
/// sequence, and that sequence itself.
struct ReadDetail {
    located: String,
    sequence: String,
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

/// Format a system time as `YYYY-MM-DD HH:MM:SS UTC`.
///
/// Done by hand rather than by pulling in a date library: the report only needs a
/// timestamp a person can read, and UTC sidesteps needing the timezone database that the
/// standard library cannot reach anyway.
fn format_utc(time: SystemTime) -> String {
    let Ok(since_epoch) = time.duration_since(SystemTime::UNIX_EPOCH) else {
        // A clock set before 1970 is not worth handling, but is not worth crashing over.
        return "unknown".to_string();
    };
    let secs = since_epoch.as_secs() as i64;
    let (days, rest) = (secs.div_euclid(86_400), secs.rem_euclid(86_400));
    let (year, month, day) = civil_from_days(days);
    let (hour, minute, second) = (rest / 3600, (rest % 3600) / 60, rest % 60);
    format!("{year:04}-{month:02}-{day:02} {hour:02}:{minute:02}:{second:02} UTC")
}

/// Days since 1970-01-01 as a calendar date, by Howard Hinnant's `civil_from_days`.
fn civil_from_days(days: i64) -> (i64, i64, i64) {
    // Shift the epoch to 0000-03-01 so that leap days fall at the end of the cycle.
    let shifted = days + 719_468;
    let era = shifted.div_euclid(146_097);
    let day_of_era = shifted.rem_euclid(146_097);
    let year_of_era =
        (day_of_era - day_of_era / 1460 + day_of_era / 36_524 - day_of_era / 146_096) / 365;
    let day_of_year = day_of_era - (365 * year_of_era + year_of_era / 4 - year_of_era / 100);
    let month_index = (5 * day_of_year + 2) / 153;
    let day = day_of_year - (153 * month_index + 2) / 5 + 1;
    let month = if month_index < 10 {
        month_index + 3
    } else {
        month_index - 9
    };
    let year = year_of_era + era * 400 + i64::from(month <= 2);
    (year, month, day)
}

/// A byte count at human scale, e.g. `1.4 MiB`.
fn human_bytes(bytes: u64) -> String {
    const UNITS: [&str; 5] = ["B", "KiB", "MiB", "GiB", "TiB"];
    let mut size = bytes as f64;
    let mut unit = 0;
    while size >= 1024.0 && unit + 1 < UNITS.len() {
        size /= 1024.0;
        unit += 1;
    }
    // Whole bytes are always whole; anything scaled reads better to one decimal place.
    if unit == 0 {
        format!("{bytes} B")
    } else {
        format!("{size:.1} {}", UNITS[unit])
    }
}

/// A count of bases at human scale, e.g. `7.6 Mbase`. Decimal rather than binary, which
/// is how sequence lengths are always quoted.
fn human_bases(bases: f64) -> String {
    const UNITS: [&str; 4] = ["base", "kbase", "Mbase", "Gbase"];
    let mut size = bases;
    let mut unit = 0;
    while size >= 1000.0 && unit + 1 < UNITS.len() {
        size /= 1000.0;
        unit += 1;
    }
    format!("{size:.1} {}", UNITS[unit])
}

/// A duration at human scale: seconds while that stays readable, then minutes and hours.
fn human_duration(elapsed: Duration) -> String {
    let secs = elapsed.as_secs_f64();
    if secs < 60.0 {
        return format!("{secs:.2} s");
    }
    let whole = elapsed.as_secs();
    let (hours, minutes, seconds) = (whole / 3600, (whole % 3600) / 60, whole % 60);
    if hours > 0 {
        format!("{hours} h {minutes:02} m {seconds:02} s")
    } else {
        format!("{minutes} m {seconds:02} s")
    }
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

    /// The read accounting as report rows. Categories that cannot arise for the given
    /// settings - the anchor ones without --start-end-segs, for instance - stay at zero
    /// and are left out, so the report only shows what actually happened.
    fn rows(&self) -> Vec<Row> {
        // Blank unless there are reads to be a share of, which also drops the column.
        let share = |n: usize| -> String {
            if self.reads() == 0 {
                String::new()
            } else {
                format!("({:.1}%)", 100.0 * n as f64 / self.reads() as f64)
            }
        };
        // The breakdown rows are shares of the reads that were not classified rather
        // than of every read, so quoting them against the same total would mislead.
        let mut rows = vec![
            Row::count(2, "reads read:", self.reads(), String::new()),
            Row::count(2, "classified:", self.classified, share(self.classified)),
            Row::count(
                2,
                "not classified:",
                self.not_classified(),
                share(self.not_classified()),
            ),
        ];
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
                rows.push(Row::count(4, label, count, String::new()));
            }
        }
        // Records that never became reads at all, so outside the accounting above.
        if self.unreadable > 0 {
            rows.push(Row::count(
                2,
                "records skipped while reading the input:",
                self.unreadable,
                String::new(),
            ));
        }
        rows
    }

    /// The read accounting on its own, laid out as [`render_rows`] describes. Production
    /// renders it as one section of [`RunReport`]; this is how the tests get at the block
    /// without building a whole report around it.
    #[cfg(test)]
    fn render(&self) -> String {
        let mut rows = vec![Row::heading("Summary")];
        rows.extend(self.rows());
        render_rows(&rows)
    }
}

/// Tally of how many reads produced each distinct segment string, written as CSV when
/// `--counts` is given.
///
/// Only classified reads are counted. Reads dropped for their length or for a missing
/// anchor never reach a classification at all and are accounted for in the run summary
/// on stderr instead. A read that was classified but carried no recognisable segments is
/// counted under the empty string: that is a real result about that read, and quite
/// different from the read having been dropped.
#[derive(Default)]
struct ClassificationCounts(HashMap<String, usize>);

impl ClassificationCounts {
    /// Record one classified read.
    fn count(&mut self, segments: &str) {
        // Only the strings actually seen are stored, so a run whose reads all classify
        // differently costs one entry per distinct result and no more.
        match self.0.get_mut(segments) {
            Some(seen) => *seen += 1,
            None => {
                self.0.insert(segments.to_string(), 1);
            }
        }
    }

    /// Render as CSV, most frequent first and then alphabetically, so that two runs over
    /// the same reads always produce the same file. A run that classified nothing still
    /// gets the header row, so the file is always valid CSV.
    fn render_csv(&self) -> String {
        let mut rows: Vec<(&String, &usize)> = self.0.iter().collect();
        rows.sort_by(|a, b| b.1.cmp(a.1).then_with(|| a.0.cmp(b.0)));
        let mut csv = String::from("segments,count\n");
        for (segments, count) in rows {
            csv.push_str(&csv_field(segments));
            csv.push(',');
            csv.push_str(&count.to_string());
            csv.push('\n');
        }
        csv
    }
}

/// A CSV field, quoted only when it has to be. Segment names come from FASTA headers, so
/// a comma or a quote in one is unlikely but perfectly legal, and left as it is would
/// shift every column after it.
fn csv_field(value: &str) -> String {
    if value.contains([',', '"', '\n', '\r']) {
        format!("\"{}\"", value.replace('"', "\"\""))
    } else {
        value.to_string()
    }
}

/// One line of the end-of-run report: a heading, a labelled count to be right-aligned
/// into the numbers column, or a labelled value that simply follows its label.
enum Row {
    Heading(String),
    Count {
        indent: usize,
        label: String,
        count: usize,
        share: String,
    },
    Value {
        indent: usize,
        label: String,
        value: String,
    },
}

impl Row {
    fn heading(title: &str) -> Row {
        Row::Heading(title.to_string())
    }

    fn count(indent: usize, label: &str, count: usize, share: String) -> Row {
        Row::Count {
            indent,
            label: label.to_string(),
            count,
            share,
        }
    }

    fn value(label: &str, value: impl Into<String>) -> Row {
        Row::Value {
            indent: 2,
            label: label.to_string(),
            value: value.into(),
        }
    }

    /// How wide this row's label is, indent included, or none for a heading.
    fn label_width(&self) -> Option<usize> {
        match self {
            Row::Heading(_) => None,
            Row::Count { indent, label, .. } | Row::Value { indent, label, .. } => {
                Some(indent + label.len())
            }
        }
    }
}

/// Lay out report rows so that every count lands in one column whatever its nesting
/// depth, and every value starts in one column too.
///
/// The columns are measured from the rows actually being printed rather than fixed, so a
/// run of a hundred reads and a run of ten million both come out square. Shares are
/// right-aligned as well, which lines up the decimal point.
fn render_rows(rows: &[Row]) -> String {
    let width = |lengths: &mut dyn Iterator<Item = usize>| lengths.max().unwrap_or(0);
    let label_width = width(&mut rows.iter().filter_map(Row::label_width));
    let count_width = width(&mut rows.iter().filter_map(|row| match row {
        Row::Count { count, .. } => Some(count.to_string().len()),
        _ => None,
    }));
    let share_width = width(&mut rows.iter().filter_map(|row| match row {
        Row::Count { share, .. } => Some(share.len()),
        _ => None,
    }));

    let mut out = String::new();
    for row in rows {
        match row {
            Row::Heading(title) => out.push_str(&format!("\n{title}\n")),
            Row::Count {
                indent,
                label,
                count,
                share,
            } => {
                let label = format!("{}{label}", " ".repeat(*indent));
                out.push_str(&format!("{label:<label_width$}  {count:>count_width$}"));
                if !share.is_empty() {
                    out.push_str(&format!("  {share:>share_width$}"));
                }
                out.push('\n');
            }
            Row::Value {
                indent,
                label,
                value,
            } => {
                let label = format!("{}{label}", " ".repeat(*indent));
                out.push_str(&format!("{label:<label_width$}  {value}\n"));
            }
        }
    }
    out
}

/// Everything the end-of-run report says about the run itself: what went in, how it was
/// configured, when it ran and how fast. The read accounting comes from [`RunSummary`].
///
/// Written to stdout so it can be captured with `> run.txt` or piped into something that
/// keeps it, separately from the warnings on stderr. The results go to their own file, so
/// nothing is competing for stdout.
struct RunReport {
    started: SystemTime,
    finished: SystemTime,
    elapsed: Duration,
    segments_file: String,
    segments_loaded: usize,
    sequences_file: String,
    sequences_format: &'static str,
    sequences_bytes: u64,
    reads: usize,
    bases: u64,
    options: Vec<(String, String)>,
}

impl RunReport {
    /// The whole report: the run, its inputs and options, the read accounting, and how
    /// fast it went. Laid out by [`render_rows`], so every column lines up across the
    /// sections rather than each being square only within itself.
    fn render(&self, summary: &RunSummary) -> String {
        let mut rows = vec![
            Row::heading("Run"),
            Row::value("started", format_utc(self.started)),
            Row::value("finished", format_utc(self.finished)),
            Row::value("elapsed", human_duration(self.elapsed)),
            Row::heading("Input"),
            Row::value("segments file", self.segments_file.clone()),
            Row::count(2, "segments loaded", self.segments_loaded, String::new()),
            Row::value("sequences file", self.sequences_file.clone()),
            Row::value("sequences format", self.sequences_format),
            Row::value("sequences size", human_bytes(self.sequences_bytes)),
            Row::heading("Options"),
        ];
        rows.extend(
            self.options
                .iter()
                .map(|(name, value)| Row::value(name, value.clone())),
        );
        rows.push(Row::heading("Summary"));
        rows.extend(summary.rows());
        rows.push(Row::heading("Throughput"));
        rows.push(Row::value("bases read", human_bases(self.bases as f64)));
        // A run so short that the clock cannot separate its ends has no meaningful rate
        // to quote, and dividing by it would print an infinity.
        let seconds = self.elapsed.as_secs_f64();
        if seconds > 0.0 {
            rows.push(Row::value(
                "reads per second",
                format!("{:.0}", self.reads as f64 / seconds),
            ));
            rows.push(Row::value(
                "bases per second",
                human_bases(self.bases as f64 / seconds),
            ));
        }
        render_rows(&rows)
    }
}

/// The auto-detected format of a sequences file.
enum InputFormat {
    Fasta,
    FastaGz,
    Fastq,
    FastqGz,
    Bam,
}

impl InputFormat {
    /// What to call this format in the run report.
    fn name(&self) -> &'static str {
        match self {
            InputFormat::Fasta => "FASTA",
            InputFormat::FastaGz => "gzipped FASTA",
            InputFormat::Fastq => "FASTQ",
            InputFormat::FastqGz => "gzipped FASTQ",
            InputFormat::Bam => "unaligned BAM",
        }
    }
}

/// Inspect a file to determine whether it's FASTA (starts with '>'), FASTQ (starts with
/// '@'), either of those gzipped, or BAM.
///
/// BAM is always BGZF-compressed (i.e. gzip-wrapped), so a gzip magic number alone doesn't
/// distinguish it from a gzipped FASTQ or FASTA file - the decompressed content is peeked
/// to tell them apart (BAM starts with the "BAM\1" magic bytes, FASTQ with '@', FASTA
/// with '>').
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
        } else if peek_n >= 1 && peek[0] == b'>' {
            Ok(InputFormat::FastaGz)
        } else {
            Err(format!(
                "Could not determine the format of '{}'. It is gzip-compressed but its contents are \
                 none of FASTQ (records start with '@'), FASTA (records start with '>') or BAM.",
                filename.display()
            ))
        }
    } else if n >= 1 && magic[0] == b'@' {
        Ok(InputFormat::Fastq)
    } else if n >= 1 && magic[0] == b'>' {
        Ok(InputFormat::Fasta)
    } else {
        Err(format!(
            "Could not determine the format of '{}'. Expected FASTQ (records start with '@'), \
             FASTA (records start with '>'), either of those gzipped, or an unaligned BAM.",
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
    Fasta(fasta::Reader<std::io::BufReader<Box<dyn Read>>>),
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
    /// What `detect_format` decided when the file was opened, kept so the run report can
    /// name it without re-reading the file - which by then might not even be there.
    format: &'static str,
    filename: std::path::PathBuf,
    bytes_read: std::sync::Arc<std::sync::atomic::AtomicU64>,
    record: fastq::Record,
    fasta_record: fasta::Record,
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
        let format = detect_format(filename)?;
        let source = match format {
            InputFormat::Fasta => {
                let raw: Box<dyn Read> = Box::new(counted(&bytes_read)?);
                Source::Fasta(fasta::Reader::new(raw))
            }
            InputFormat::FastaGz => {
                let raw: Box<dyn Read> =
                    Box::new(flate2::read::MultiGzDecoder::new(counted(&bytes_read)?));
                Source::Fasta(fasta::Reader::new(raw))
            }
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
            format: format.name(),
            filename: filename.to_path_buf(),
            bytes_read,
            record: fastq::Record::new(),
            fasta_record: fasta::Record::new(),
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
            // FASTA carries no qualities, which the classifier never looks at anyway, so
            // apart from the record type this mirrors the FASTQ branch below exactly.
            Source::Fasta(reader) => loop {
                match reader.read(&mut self.fasta_record) {
                    Ok(()) => {
                        if self.fasta_record.is_empty() {
                            return Ok(None);
                        }
                        self.index += 1;
                        self.consecutive_errors = 0;
                        if self.fasta_record.seq().is_empty() {
                            warn(format!(
                                "skipping read '{}' in '{}': it has no sequence",
                                self.fasta_record.id(),
                                self.filename.display()
                            ));
                            self.skipped += 1;
                            continue;
                        }
                        return Ok(Some(RawRead {
                            name: self.fasta_record.id().to_string(),
                            seq: self.fasta_record.seq().to_ascii_uppercase(),
                        }));
                    }
                    Err(e) => {
                        self.index += 1;
                        self.skipped += 1;
                        self.consecutive_errors += 1;
                        // The name is only known if the parser reached the header line.
                        if self.fasta_record.id().is_empty() {
                            warn(format!(
                                "skipping malformed record {} in '{}': {e}",
                                self.index,
                                self.filename.display()
                            ));
                        } else {
                            warn(format!(
                                "skipping malformed read '{}' in '{}': {e}",
                                self.fasta_record.id(),
                                self.filename.display()
                            ));
                        }
                        if self.consecutive_errors >= MAX_CONSECUTIVE_BAD_RECORDS {
                            return Err(format!(
                                "Gave up on '{}' after {} malformed records in a row. The file \
                                 is probably truncated or not FASTA at all.",
                                self.filename.display(),
                                self.consecutive_errors
                            ));
                        }
                    }
                }
            },
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
/// Bases match when their IUPAC nucleotide sets intersect, so an ambiguity code in an
/// anchor scores a full match against any base it permits. See [`BASE_MASK`].
///
/// Mirrors the affine gap model of the library aligner it replaces, where a run of `k`
/// gaps costs `gap_open + gap_extend * k` - so with both set to -2 a single gap costs -4,
/// not -2. Three running rows track alignments ending in a match, in a gap in the target,
/// and in a gap in the query, each carrying the target position its path began at.
fn best_semiglobal(query: &[u8], target: &[u8], scoring: (i32, i32, i32, i32)) -> AnchorAlignment {
    let (match_score, mismatch_score, gap_open, gap_extend) = scoring;
    let query_len = query.len();
    let target_len = target.len();
    // As in `align_multiple`: one pass to encode, then set intersection per cell.
    let query_masks = base_masks(query);
    let target_masks = base_masks(target);
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
        let query_mask = query_masks[qi - 1];
        for j in 1..=target_len {
            let sub = if query_mask & target_masks[j - 1] != 0 {
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

/// The 'start' and 'end' anchors with their reverse complements and the raw score each
/// has to reach, all computed once per run. The minimums come from each anchor's own
/// threshold, so an anchor may be held to a different standard from the segments between
/// them - and from each other.
struct AnchorSeqs {
    start: Vec<u8>,
    start_rc: Vec<u8>,
    end: Vec<u8>,
    end_rc: Vec<u8>,
    start_min: i32,
    end_min: i32,
}

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

/// One segment found in a read: the name it is reported under, with a trailing `*` if it
/// was found on the reverse strand, and the half-open span of the extracted sequence it
/// covers.
struct FoundSegment {
    name: String,
    start: usize,
    end: usize,
}

/// The segments found in a read, in the order they occur along it.
///
/// With `start_end_segs` the read has already been trimmed to run from 'start' to 'end',
/// so those two are not searched for within it; [`segment_string`] adds them back as the
/// bounds. Without it the read is classified as given and a segment named 'start' or
/// 'end' is treated like any other.
///
/// With `omit_rc_segs` each segment is looked for on the forward strand only, so no hit
/// can ever be reported under a starred name. That is the caller asserting the strand
/// every segment sits on is already known - which anchoring establishes, since the read
/// has been reoriented by then - and it halves the alignment work. Without it both
/// strands are searched, which is the default and what to keep when the orientation of a
/// segment within the construct is not known in advance.
///
/// Each segment is found at its own threshold, resolved when the segments file was
/// loaded, so there is no global score to apply here.
fn classify_read_segments(
    segments: &HashMap<String, Segment>,
    read_seq: &[u8],
    start_end_segs: bool,
    omit_rc_segs: bool,
) -> Vec<FoundSegment> {
    let aligner = MultiTracebackAligner::new(
        2,  // match score
        -1, // mismatch score
        -2, // gap score
    );
    let mut all_alignments: Vec<MultiAlignment> = Vec::new();
    for (seg_name, segment) in segments {
        let seg_seq = &segment.seq;
        if start_end_segs && (seg_name == "start" || seg_name == "end") {
            continue;
        }
        // Rejected when the segments file is loaded, but guard anyway: a zero-length
        // segment scores zero everywhere, so it would "match" at every position.
        if seg_seq.is_empty() {
            continue;
        }
        let min_score = (seg_seq.len() as f32 * segment.min_norm_score) as i32;
        all_alignments.append(&mut aligner.align_multiple(seg_name, seg_seq, read_seq, min_score));
        if !omit_rc_segs {
            all_alignments.append(&mut aligner.align_multiple(
                &format!("{}*", seg_name),
                &dna::revcomp(seg_seq),
                read_seq,
                min_score,
            ));
        }
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
    all_alignments
        .into_iter()
        .filter(|aln| !aln.filter)
        .map(|aln| FoundSegment {
            name: aln.name,
            start: aln.target_start,
            end: aln.target_end,
        })
        .collect()
}

/// The segment string: the names in read order joined by `-`.
fn segment_string(found: &[FoundSegment]) -> String {
    found
        .iter()
        .map(|segment| segment.name.as_str())
        .collect::<Vec<_>>()
        .join("-")
}

/// The same segments with the span each covers, as `seg1[1:13],seg2*[15:20]`.
///
/// Entries are separated by commas and the two positions within one by a colon, so the
/// column can be split on `,` without the separators being confused for one another.
///
/// Positions are 1-based and inclusive of both ends, counted along the extracted sequence
/// written beside them rather than along the read as it arrived - with `--start-end-segs`
/// the two differ, since the read has been trimmed and possibly reoriented by then.
fn located_segment_string(found: &[FoundSegment]) -> String {
    found
        .iter()
        .map(|segment| format!("{}[{}:{}]", segment.name, segment.start + 1, segment.end))
        .collect::<Vec<_>>()
        .join(",")
}

/// Classify one extracted sequence and package it for output.
///
/// `clean_seq` is the read as it is actually classified: the whole read when no anchors
/// were given, and the trimmed, reoriented span between them when they were. That is the
/// sequence `--detailed-output` writes and the one the positions beside it are counted
/// along, so the three columns always describe the same thing.
///
/// `anchor_spans` carries how much of `clean_seq` the 'start' and 'end' anchors take up,
/// and is `None` when the read was classified without anchors. The extracted sequence
/// runs from the first base of 'start' to the last base of 'end' whichever strand the
/// read arrived on, so the anchors are always its two ends and their spans follow from
/// their lengths alone.
fn classified_read(
    read_name: String,
    clean_seq: &[u8],
    segments: &HashMap<String, Segment>,
    anchor_spans: Option<(usize, usize)>,
    omit_rc_segs: bool,
    detailed: bool,
) -> SegmentedRead {
    let mut found =
        classify_read_segments(segments, clean_seq, anchor_spans.is_some(), omit_rc_segs);
    if let Some((start_span, end_span)) = anchor_spans {
        // Put back at their known positions rather than searched for: they were located
        // once already, to trim the read to them.
        found.insert(
            0,
            FoundSegment {
                name: "start".to_string(),
                start: 0,
                end: start_span,
            },
        );
        found.push(FoundSegment {
            name: "end".to_string(),
            // Saturating because a pathological read could align both anchors across
            // more bases than lie between them; better a clamped span than a panic.
            start: clean_seq.len().saturating_sub(end_span),
            end: clean_seq.len(),
        });
    }
    SegmentedRead {
        name: read_name,
        segments: segment_string(&found),
        detail: detailed.then(|| ReadDetail {
            located: located_segment_string(&found),
            sequence: String::from_utf8_lossy(clean_seq).into_owned(),
        }),
    }
}

/// Classify segments across all reads loaded from a sequences file (FASTQ or BAM).
fn process_reads(
    records: Vec<RawRead>,
    segments: &HashMap<String, Segment>,
    // Minimum and maximum read length, where zero on either side means no bound there.
    len_check: (usize, usize),
    circular: bool,
    start_end_segs: bool,
    omit_rc_segs: bool,
    detailed: bool,
) -> Result<(Vec<SegmentedRead>, RunSummary)> {
    // Precompute start/end revcomps and score minimums once rather than per-read
    let start_end_seqs: Option<AnchorSeqs> = if start_end_segs {
        let anchor = |name: &str| -> Result<Segment> {
            segments.get(name).cloned().ok_or_else(|| {
                format!(
                    "--start-end-segs was given, but the segments file has no segment named \
                     '{name}'. Add it, or drop --start-end-segs to classify reads as they are."
                )
            })
        };
        let start = anchor("start")?;
        let end = anchor("end")?;
        Some(AnchorSeqs {
            start_rc: dna::revcomp(&start.seq),
            end_rc: dna::revcomp(&end.seq),
            start_min: (start.seq.len() as f32 * start.min_norm_score) as i32,
            end_min: (end.seq.len() as f32 * end.min_norm_score) as i32,
            start: start.seq,
            end: end.seq,
        })
    } else {
        None
    };

    let clean_seqs: Vec<Outcome> = records
        .par_iter()
        .map(|record| {
            let read_name = record.name.clone();
            let read_seq = record.seq.clone();

            let read_len = read_seq.len();
            // Zero means no bound on that side, independently of the other, so a
            // minimum can be set without also having to invent a maximum.
            if len_check.0 > 0 && read_len < len_check.0 {
                return Outcome::Rejected(Rejected::TooShort);
            }
            if len_check.1 > 0 && read_len > len_check.1 {
                return Outcome::Rejected(Rejected::TooLong);
            }

            if !start_end_segs {
                // No anchors: classify the read as it arrived, on whichever strand.
                let clean_seq = cut_reorient_seq(&read_seq, 0, read_seq.len(), false);
                Outcome::Classified(classified_read(
                    read_name,
                    &clean_seq,
                    segments,
                    None,
                    omit_rc_segs,
                    detailed,
                ))
            } else if let Some(anchors) = start_end_seqs.as_ref() {
                let anchor = |seq: &[u8]| best_semiglobal(seq, &read_seq, ANCHOR_SCORING);
                let align_start = anchor(&anchors.start);
                let align_end = anchor(&anchors.end);
                let align_start_rc = anchor(&anchors.start_rc);
                let align_end_rc = anchor(&anchors.end_rc);

                let (start_min, end_min) = (anchors.start_min, anchors.end_min);
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
                        Outcome::Classified(classified_read(
                            read_name,
                            &clean_seq,
                            segments,
                            Some((
                                align_start.yend - align_start.ystart,
                                align_end.yend - align_end.ystart,
                            )),
                            omit_rc_segs,
                            detailed,
                        ))
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
                        Outcome::Classified(classified_read(
                            read_name,
                            &clean_seq,
                            segments,
                            Some((
                                align_start_rc.yend - align_start_rc.ystart,
                                align_end_rc.yend - align_end_rc.ystart,
                            )),
                            omit_rc_segs,
                            detailed,
                        ))
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
/// the best one. Matching is by IUPAC nucleotide set, so segments may be degenerate.
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
    /// Bases match when their IUPAC nucleotide sets intersect, so an ambiguity code in a
    /// segment scores a full match against any base it permits. See [`BASE_MASK`].
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
        // Encoded once here rather than per cell: the encoding is linear in the sequence
        // while the matrix below is quadratic, so it disappears into the noise.
        let query_masks = base_masks(query);
        let target_masks = base_masks(target);

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
            let query_mask = query_masks[i - 1];
            for j in 1..=target_len {
                let match_score = if query_mask & target_masks[j - 1] != 0 {
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

    /// Minimum read length to process. Zero means no minimum, which is the default: no
    /// read is dropped for its length unless a bound is asked for
    #[arg(long, default_value_t = 0)]
    min_seq_len: usize,

    /// Maximum read length to process. Zero means no maximum, which is the default
    #[arg(long, default_value_t = 0)]
    max_seq_len: usize,

    /// Use the 'start' and 'end' segments to anchor, trim and orient each read. They
    /// then bound the segment string rather than being reported within it. Without
    /// this, reads are classified as given and any segment named 'start' or 'end' is
    /// treated as an ordinary segment
    #[arg(short = 'd', long, default_value_t = false)]
    start_end_segs: bool,

    /// Minimum normalised score, between 0.0 and 2.0 (perfect match). Applies to every
    /// segment that does not name its own with --per-segment-scores
    #[arg(short = 's', long, default_value_t = 1.5)]
    min_norm_score: f32,

    /// Read a per-segment minimum normalised score from square brackets after each
    /// segment's name, so that '>SEG1[1.7]' is the segment 'SEG1' found at 1.7. Segments
    /// without brackets use --min-norm-score. Without this flag the brackets are part of
    /// the name, as they were before the option existed
    #[arg(long, default_value_t = false)]
    per_segment_scores: bool,

    /// Look for each segment on the forward strand only, so nothing is reported under a
    /// starred name. Use it when every segment's orientation within the construct is
    /// already known - anchoring with --start-end-segs establishes it - to halve the
    /// alignment work. Read orientation is unaffected: anchors are still matched on both
    /// strands and reversed reads still recovered
    #[arg(long, default_value_t = false)]
    omit_rc_segs: bool,

    /// Sequences to classify are circular (i.e., from plasmids)
    #[arg(short = 'c', long, default_value_t = false)]
    circular: bool,

    /// Threads to use
    #[arg(short = 't', long, default_value_t = 1)]
    threads: usize,

    /// Where to write the per-read classifications
    #[arg(short = 'o', long, default_value = "classifications.txt")]
    classifications: String,

    /// Also write a CSV of every distinct classification and how many reads produced it,
    /// most frequent first
    #[arg(long)]
    counts: Option<String>,

    /// Add two columns to the classifications file: where each segment sits in the
    /// extracted sequence, as 'seg1[1,13],seg2*[15,20]', and that sequence itself
    #[arg(long, default_value_t = false)]
    detailed_output: bool,
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
    // Both clocks: the wall clock for the timestamps a person reads, and a monotonic one
    // for the elapsed time, which must not be thrown off by the clock being adjusted.
    let started = SystemTime::now();
    let start_instant = Instant::now();

    if !(0.0..=2.0).contains(&cli.min_norm_score) {
        return Err(format!(
            "--min-norm-score must be between 0.0 and 2.0 (2.0 is a perfect match), but {} was given.",
            cli.min_norm_score
        ));
    }
    // Only meaningful when there is an upper bound: with --max-seq-len at its default of
    // zero there is nothing for the minimum to exceed.
    if cli.max_seq_len > 0 && cli.min_seq_len > cli.max_seq_len {
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

    let segments = load_fasta(
        Path::new(&cli.segments),
        cli.min_norm_score,
        cli.per_segment_scores,
    )?;
    for (name, chance, required) in degenerate_segments(&segments) {
        warn(format!(
            "segment '{name}' is ambiguous enough to match {:.0}% of random bases, which is at or \
             above the {:.0}% identity its minimum normalised score asks for, so expect to find it \
             almost anywhere. Use a more specific sequence, or raise the score it is found at.",
            chance * 100.0,
            required * 100.0,
        ));
    }
    let sequences = Path::new(&cli.sequences);
    let mut stream = ReadStream::open(sequences)?;

    let f = File::create(Path::new(&cli.classifications)).map_err(|e| {
        format!(
            "Could not create classifications file '{}': {e}",
            cli.classifications
        )
    })?;
    let mut f = BufWriter::new(f);

    // Created up front alongside the classifications file rather than at the end of the
    // run, so an unwritable path is reported before the work is done instead of after it.
    let mut counts_file = match &cli.counts {
        Some(path) => Some((
            path.clone(),
            File::create(Path::new(path))
                .map_err(|e| format!("Could not create counts file '{path}': {e}"))?,
        )),
        None => None,
    };

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
    let mut classifications = ClassificationCounts::default();
    let mut reads_seen = 0usize;
    let mut bases_seen = 0u64;
    loop {
        let chunk = stream.next_chunk(MAX_CHUNK_BASES, MAX_CHUNK_READS)?;
        if chunk.is_empty() {
            break;
        }
        reads_seen += chunk.len();
        bases_seen += chunk
            .iter()
            .map(|record| record.seq.len() as u64)
            .sum::<u64>();
        let (clean_seqs, chunk_summary) = process_reads(
            chunk,
            &segments,
            (cli.min_seq_len, cli.max_seq_len),
            cli.circular,
            cli.start_end_segs,
            cli.omit_rc_segs,
            cli.detailed_output,
        )?;
        for r in clean_seqs {
            let row = match &r.detail {
                Some(detail) => format!(
                    "{}\t{}\t{}\t{}",
                    r.name, r.segments, detail.located, detail.sequence
                ),
                None => format!("{}\t{}", r.name, r.segments),
            };
            writeln!(f, "{row}").map_err(|e| {
                format!(
                    "Could not write to classifications file '{}': {e}",
                    cli.classifications
                )
            })?;
            // Counted only when asked for: on a run where every read classifies
            // differently the tally grows with the input, and there is no reason to pay
            // for that when no tally was requested.
            if counts_file.is_some() {
                classifications.count(&r.segments);
            }
        }
        summary.merge(chunk_summary);
        // Advanced once the chunk is written, so the bar means "done", not "read ahead".
        progress.set_position(stream.bytes_read().min(total_bytes));
        progress.set_message(format!("{reads_seen} reads"));
    }
    progress.set_position(total_bytes);
    progress.finish_and_clear();

    f.flush().map_err(|e| {
        format!(
            "Could not finish writing classifications file '{}': {e}",
            cli.classifications
        )
    })?;

    if let Some((path, file)) = counts_file.as_mut() {
        file.write_all(classifications.render_csv().as_bytes())
            .map_err(|e| format!("Could not write counts file '{path}': {e}"))?;
    }

    report_skipped(stream.skipped, reads_seen, sequences);
    if reads_seen == 0 {
        warn(format!(
            "no usable reads were loaded from '{}'",
            cli.sequences
        ));
    }
    summary.unreadable = stream.skipped;

    let on_off = |enabled: bool| if enabled { "on" } else { "off" }.to_string();
    let report = RunReport {
        started,
        finished: SystemTime::now(),
        elapsed: start_instant.elapsed(),
        segments_file: cli.segments.clone(),
        segments_loaded: segments.len(),
        sequences_file: cli.sequences.clone(),
        sequences_format: stream.format,
        sequences_bytes: total_bytes,
        reads: reads_seen,
        bases: bases_seen,
        // Every option at the value it actually ran with, so the report is a record of
        // what produced the results beside it and not just of how they turned out.
        options: vec![
            (
                "--min-norm-score".to_string(),
                cli.min_norm_score.to_string(),
            ),
            (
                "--per-segment-scores".to_string(),
                on_off(cli.per_segment_scores),
            ),
            ("--start-end-segs".to_string(), on_off(cli.start_end_segs)),
            ("--omit-rc-segs".to_string(), on_off(cli.omit_rc_segs)),
            ("--circular".to_string(), on_off(cli.circular)),
            ("--min-seq-len".to_string(), cli.min_seq_len.to_string()),
            ("--max-seq-len".to_string(), cli.max_seq_len.to_string()),
            ("--threads".to_string(), cli.threads.to_string()),
            ("--detailed-output".to_string(), on_off(cli.detailed_output)),
            ("--classifications".to_string(), cli.classifications.clone()),
            (
                "--counts".to_string(),
                cli.counts
                    .clone()
                    .unwrap_or_else(|| "(not written)".to_string()),
            ),
        ],
    };
    print!("{}", report.render(&summary));
    Ok(())
}
