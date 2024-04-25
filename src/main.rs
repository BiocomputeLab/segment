// Segment command line tool
// Author: Thomas Gorochowski <tom@chofski.co.uk>

use bio;
use bio::alignment::pairwise::*;
use bio::alignment::Alignment;
use bio::alphabets::dna;
use bio::io::fasta;
use bio::io::fastq;

use std::collections::HashMap;
use std::path::Path;

/// Loads a FASTA file containing the `start` and `end` sequences plus and
/// other segment sequences that should be used when classifying a read.
fn load_fasta(filename: &Path) -> HashMap<String, Vec<u8>> {
    let mut fasta_data: HashMap<String, Vec<u8>> = HashMap::new();
    let reader = fasta::Reader::from_file(filename).unwrap();
    for result in reader.records() {
        let record = result.expect("Error during FASTA record parsing");
        fasta_data.insert(record.id().to_string(), record.seq().to_owned());
    }
    fasta_data
}

struct SegmentedRead {
    name: String,
    seq: Vec<u8>,
    qual: Vec<u8>,
    segments: String,
}

fn cut_reorient_seq(seq: &Vec<u8>, start_idx: usize, end_idx: usize, revcomp: bool) -> Vec<u8> {
    let new_seq: Vec<u8>;
    if start_idx < end_idx {
        new_seq = seq[start_idx..end_idx].to_vec();
    } else {
        if revcomp == true {
            let rc_seq = dna::revcomp(seq);
            new_seq = rc_seq[(seq.len() - start_idx)..(seq.len() - end_idx)].to_vec();
        } else {
            let mut rev_seq = seq.clone();
            rev_seq.reverse();
            new_seq = rev_seq[(seq.len() - start_idx)..(seq.len() - end_idx)].to_vec();
        }
    }
    new_seq
}

fn align_segment(segment_seq: &Vec<u8>, read_seq: &Vec<u8>) -> Alignment {
    let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
    // gap open score: -5, gap extension score: -1
    let mut aligner = Aligner::with_capacity(segment_seq.len(), read_seq.len(), -5, -1, &score);
    let alignment = aligner.semiglobal(segment_seq, read_seq);
    alignment
}

/// Calculate a reasonable alignment mapping score to trust that a segment has been found.
fn score_threshold(read_len: usize) -> i32 {
    let read_len_i32 = read_len as i32;
    let min_score = read_len_i32 - (5 * (read_len_i32 / 14));
    println!("{}", min_score);
    min_score
}

fn classify_read_segments(segments: &HashMap<String, Vec<u8>>, read_seq: &Vec<u8>) -> String {
    let start_str = "start";
    let mut seg_str = start_str.to_string();
    for (seg_name, seg_seq) in segments {
        if seg_name == "start" || seg_name == "end" {
            continue;
        }
        let seg_alignment = align_segment(seg_seq, read_seq);
        let min_score = score_threshold(seg_seq.len());
        if seg_alignment.score >= min_score {
            seg_str.push_str(" - ");
            seg_str.push_str(seg_name.as_str());
        } else {
            let seg_seq_rc = dna::revcomp(seg_seq).clone();
            let seg_alignment_rc = align_segment(&seg_seq_rc, read_seq);
            if seg_alignment_rc.score >= min_score {
                let seg_name_rc: String = format!("{}*", seg_name);
                seg_str.push_str(" - ");
                seg_str.push_str(seg_name_rc.as_str());
            }
        }
    }
    seg_str.push_str(" - end");
    seg_str
}

fn process_fastq(
    filename: &Path,
    segments: &HashMap<String, Vec<u8>>,
    len_check: (usize, usize),
) -> Vec<SegmentedRead> {
    let mut clean_seqs: Vec<SegmentedRead> = Vec::new();
    let reader = fastq::Reader::from_file(filename).unwrap();
    let mut read_num = 0;
    for result in reader.records() {
        if read_num % 100 == 0 {
            println!("Processing read: {}", read_num);
        }
        read_num += 1;

        let record = result.expect("Error during FASTQ record parsing");
        let read_name = record.id().to_string();
        let read_seq = record.seq().to_owned();
        let read_qual = record.qual().to_owned();

        // Check the length of the read is in expected bounds.
        let read_len = read_seq.len();
        if len_check != (0, 0) && (read_len < len_check.0 || read_len > len_check.1) {
            continue;
        }

        // Find the start and end segment in the read (if it isn't found then don't process the
        // read further).
        let start_seq = match segments.get("start") {
            Some(value) => value,
            None => panic!("No 'start' segment sequence found."),
        };
        let start_rc_seq = dna::revcomp(start_seq);
        let end_seq = match segments.get("end") {
            Some(value) => value,
            None => panic!("No 'end' segment sequence found."),
        };
        let end_rc_seq = dna::revcomp(end_seq);
        // Perform the alignments
        // TODO: could reduce this to stop when threshold score met
        let align_start = align_segment(start_seq, &read_seq);
        let align_start_rc = align_segment(&start_rc_seq, &read_seq);
        let align_end = align_segment(end_seq, &read_seq);
        let align_end_rc = align_segment(&end_rc_seq, &read_seq);

        if align_start.score > align_start_rc.score && align_end.score > align_end_rc.score {
            // Read in the correct orientation
            if align_start.score >= score_threshold(start_seq.len())
                && align_end.score >= score_threshold(end_seq.len())
            {
                let clean_seq =
                    cut_reorient_seq(&read_seq, align_start.ystart, align_end.yend, true);
                let clean_qual =
                    cut_reorient_seq(&read_qual, align_start.ystart, align_end.yend, false);
                let seg_str = classify_read_segments(segments, &clean_seq);
                let seg_data = SegmentedRead {
                    name: read_name,
                    seq: clean_seq,
                    qual: clean_qual,
                    segments: seg_str,
                };
                clean_seqs.push(seg_data);
            }
        } else {
            // Read in the reverse orientation
            if align_start_rc.score >= score_threshold(start_seq.len())
                && align_end_rc.score >= score_threshold(end_seq.len())
            {
                let clean_seq =
                    cut_reorient_seq(&read_seq, align_start_rc.yend, align_end_rc.ystart, true);
                let clean_qual =
                    cut_reorient_seq(&read_qual, align_start_rc.yend, align_end_rc.ystart, false);
                let seg_str = classify_read_segments(segments, &clean_seq);
                let seg_data = SegmentedRead {
                    name: read_name,
                    seq: clean_seq,
                    qual: clean_qual,
                    segments: seg_str,
                };
                clean_seqs.push(seg_data);
            }
        }
    }
    clean_seqs
}

fn main() {
    let fasta_path = Path::new("./data/refs.fasta");
    let fastq_path = Path::new("./data/reads.fastq");

    let segments = load_fasta(fasta_path);
    for (name, seq) in &segments {
        let s = String::from_utf8(seq.clone()).unwrap();
        println!("Name: {}\nSequnce: {}", name, s);
    }

    let clean_seqs = process_fastq(&fastq_path, &segments, (900, 1200));
    println!("{:?}", clean_seqs.len());
    for r in clean_seqs {
        let s = String::from_utf8(r.seq.clone()).unwrap();
        println!("Name: {}\nSequnce: {}\nSegments: {}", r.name, s, r.segments);
    }
}
