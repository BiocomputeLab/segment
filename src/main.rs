// Segment command line tool
// Author: Thomas E. Gorochowski <tom@chofski.co.uk>

// We have some items that are not used, don't let them raise warnings
// when we are compiling.
#![allow(unused_variables, dead_code)]

use std::env;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;
use std::collections::HashMap;
use bio;
use bio::alignment::{pairwise::*, Alignment};
use bio::alphabets::dna;
use bio::io::{fasta, fastq};

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

/// Structure to hold the key sequence and segment information about a read.
struct SegmentedRead {
    name: String,
    seq: Vec<u8>,
    qual: Vec<u8>,
    segments: String,
}

/// Cut out the relevant region of a read and generate the reverse complement
/// so that a sequence from 'start' to 'end' segments is returned.
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

/// Carry out a semi-global alignment of the segment against a read.
fn align_segment(segment_seq: &Vec<u8>, read_seq: &Vec<u8>) -> Alignment {
    let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
    // gap open score: -5, gap extension score: -1
    let mut aligner = Aligner::with_capacity(segment_seq.len(), read_seq.len(), -1, -1, &score);
    let alignment = aligner.semiglobal(segment_seq, read_seq);
    alignment
}

/// Calculate a reasonable alignment mapping score to trust that a segment has been found.
fn score_threshold(read_len: usize) -> i32 {
    let read_len_i32 = read_len as i32;
    // Threshold is that there is at least a 60% match to the reference
    let min_score = (read_len_i32 as f32 * 0.6) as i32;
    min_score
}

/// Classify which segments and their orientations are present within a read.
fn classify_read_segments(segments: &HashMap<String, Vec<u8>>, read_seq: &Vec<u8>) -> String {
    let start_str = "start";
    let mut seg_str = start_str.to_string();
    let mut segs: Vec<(usize, String)> = Vec::new();
    for (seg_name, seg_seq) in segments {
        if seg_name == "start" || seg_name == "end" {
            continue;
        }
        let seg_alignment = align_segment(seg_seq, read_seq);
        let min_score = score_threshold(seg_seq.len());
        if seg_alignment.score >= min_score {
            segs.push((seg_alignment.ystart, seg_name.clone()))
        } else {
            let seg_seq_rc = dna::revcomp(seg_seq).clone();
            let seg_alignment_rc = align_segment(&seg_seq_rc, read_seq);
            if seg_alignment_rc.score >= min_score {
                let seg_name_rc: String = format!("{}*", seg_name);
                segs.push((seg_alignment_rc.ystart, seg_name_rc))
            }
        }
    }
    // Sort the segments found on start position
    segs.sort_by_key(|tuple| tuple.0);
    // Generate the segment string
    for (seg_pos, seg_name) in segs {
        seg_str.push_str("-");
        seg_str.push_str(seg_name.as_str());
    }
    seg_str.push_str("-end");
    seg_str
}

/// Load the data from a FASTQ file and extract sequences from each read
/// corresponding to the presence of both a 'start' and 'end' segment.
fn process_fastq(
    filename: &Path,
    segments: &HashMap<String, Vec<u8>>,
    len_check: (usize, usize),
) -> Vec<SegmentedRead> {
    let mut clean_seqs: Vec<SegmentedRead> = Vec::new();
    let reader = fastq::Reader::from_file(filename).unwrap();
    let mut read_num = 0;
    for result in reader.records() {
        if read_num % 1000 == 0 {
            println!("Processing read: {}", read_num);
        }
        read_num += 1;
        // Gather information about the read
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
        let align_start = align_segment(&start_seq, &read_seq);
        let align_end = align_segment(&end_seq, &read_seq);
        let align_start_rc = align_segment(&start_rc_seq, &read_seq);
        let align_end_rc = align_segment(&end_rc_seq, &read_seq);
        // Check to see if the read is valid (contains start and end) and reorient if needed
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

/// Main function to process command line arguments and perform the analysis.
/// TODO: Update to use clap (https://docs.rs/clap/latest/clap/)
fn main() {
    let args: Vec<String> = env::args().collect();
    let fasta_path = Path::new(&args[1]);
    let fastq_path = Path::new(&args[2]);
    let min_read_len = args[3].parse::<usize>().unwrap();
    let max_read_len = args[4].parse::<usize>().unwrap();
    let output_path = Path::new(&args[5]);
    let segments = load_fasta(fasta_path);
    let clean_seqs = process_fastq(&fastq_path, &segments, (min_read_len, max_read_len));
    // Write the classified segments to file
    let f = File::create(output_path).expect("unable to create file");
    let mut f = BufWriter::new(f);
    for r in clean_seqs {
        //let s = String::from_utf8(r.seq.clone()).unwrap();
        writeln!(f, "{}\t{}", r.name, r.segments).expect("unable to write");
    }
    f.flush().unwrap();
}
