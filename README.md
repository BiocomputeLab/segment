# segment

A flexible analysis tool for classifying the order and orientation of sequence segments in genetic parts, circuits and libraries.

Given a FASTA file of segment sequences (e.g., attachment sites, barcodes, promoters, parts – anything you can write down) and a file of long sequencing reads, `segment` aligns every segment against every read and reports, for each read, which segments it contains, in what order, and in what orientation. It was written for nanopore reads of recombinase-based constructs, but nothing about it is specific to that application.

## Installation

To use `segment` you will need to first compile the source code contained in this repository. This requires that the `rustc` compiler is available on your system (version 1.85 or later, for the 2024 edition). Further details about how to install `rustc` can be found at the following URL and we recommend the use of `rustup` to maintain your environment and ensure all other supporting tools are available.

https://www.rust-lang.org/tools/install

Once `rustc` is installed, compilation of `segment` is performed by running the following command from within the root of the code repository:

```sh
cargo build --release
```

This should take a couple minutes to complete and the generated `segment` executable can be found in the `target/release` directory. We recommend placing the `segment` executable in a location found in your `PATH` environment variable to make running the tool easier from the command line. Alternatively, `cargo` can do this for you:

```sh
cargo install --path .
```

It should be possible to compile and run `segment` on any system for which `rustc` is available (e.g., Windows, MacOS and Linux). All dependencies are declared in `Cargo.toml` and are fetched automatically by `cargo`. Testing of the software has been performed exclusively on MacOS Sequoia (15.5).

## Usage

To run `segment` it has the following usage:

```sh
segment --segments REF_FASTA --sequences READS [OPTIONS] --classifications OUTPUT_FILE
```

where `REF_FASTA` is a file containing the reference segments in FASTA format, `READS` is the sequencing data (FASTA, FASTQ, either of those gzipped, or unaligned BAM), and `OUTPUT_FILE` is the filename to write the results to.

By default every read is classified exactly as it arrives. Passing `-d/--start-end-segs` instead uses the segments named `start` and `end` to anchor each read: reads are oriented, trimmed to the region between the two anchors, and dropped if either anchor is missing. This is usually what you want when the reads contain adapters or other flanking sequence.

### Required arguments

| Argument | Description |
|---|---|
| `--segments <FILE>` | FASTA file of segment sequences to search for, one record per segment. Record IDs become the names used in the output. Names must be unique and no record may be empty. |
| `--sequences <FILE>` | Reads to classify: FASTA, FASTQ, either of those gzipped, or unaligned BAM. The format is detected from the file contents, so the extension is irrelevant. |

### Optional arguments

| Argument | Default | Description |
|---|---|---|
| `-o`, `--classifications <FILE>` | `classifications.txt` | Where to write the tab-separated per-read results. Overwritten if it exists. |
| `--counts <FILE>` | off | Also write a CSV counting how many reads produced each distinct classification. See [Classification counts](#classification-counts). |
| `--detailed-output` | off | Add two columns to the classifications file: where each segment sits in the extracted sequence, and that sequence itself. See [Detailed output](#detailed-output). |
| `-d`, `--start-end-segs` | off | Use the segments named `start` and `end` to anchor, trim and orient each read; they then bound the segment string instead of being reported within it. Reads missing either anchor are dropped. Errors if the segments file has no `start` or `end` record. |
| `-c`, `--circular` | off | Treat reads as coming from a circular template, so a read that crosses the origin (`end` before `start`) is rotated back into order rather than rejected. Only meaningful together with `-d`. |
| `-s`, `--min-norm-score <FLOAT>` | `1.5` | Minimum normalised alignment score for a segment to count as found, from `0.0` to `2.0`, where `2.0` is a perfect match. See [Choosing a score threshold](#choosing-a-score-threshold). |
| `--per-segment-scores` | off | Use each segment's own threshold from square brackets after its name, as in `>SEG1[1.7]`. Segments without brackets use `--min-norm-score`. The brackets are stripped from the name either way. See [Per-segment thresholds](#per-segment-thresholds). |
| `--min-seq-len <INT>` | `0` | Reads shorter than this are dropped before alignment. `0` means no minimum, so by default no read is dropped for being short. |
| `--max-seq-len <INT>` | `0` | Reads longer than this are dropped before alignment. `0` means no maximum. When both are set, this must not be less than `--min-seq-len`. |
| `-t`, `--threads <INT>` | `1` | Reads are classified in parallel across this many threads. |
| `-h`, `--help` | | Print help. |
| `-V`, `--version` | | Print version. |

Peak memory does not grow with the size of the input. Reads are streamed through in bounded chunks, and all alignment runs in space linear in the read length and independent of segment length, so long segments and long anchors cost no more than short ones. The practical ceiling is a single very long read, not the file. `--detailed-output` roughly doubles what a chunk holds, since it carries each read's extracted sequence until the chunk is written, but the chunk is still bounded.

## Example

To run `segment` on the small example data set provided in the `data` directory, we would therefore use the command (running from the root of the code repository):

```sh
segment --segments data/refs.fasta --sequences data/reads.fastq \
        --start-end-segs --min-seq-len 500 --max-seq-len 2000 \
        --classifications classifications.txt
```

This should generate a file classifying each read in the input FASTQ data file that passes the read length filtering and in which both anchors are found. The first few lines are:

```
ad6f445e-17e6-49ce-bf40-67da322006c8	start-attP_Bxb1*-end
27bcf6c4-9226-4fa9-af0a-1ed4294d8768	start-end
715b0122-721e-4dc1-acff-ca723162392c	start-attP_Bxb1*-end
fe1b90af-4e41-4c33-ab87-eee17512462c	start-attP_Bxb1*-end
cb24e39a-8ee8-4988-a788-250ccd99349c	start-attB_Tp901-attB_Int5-end
c689dd31-4e98-4154-ab05-a8b0abed9436	start-end
8c19bcd4-b228-4e0d-a257-e5116a0a4322	start-attP_Bxb1*-end
732ad32d-6079-4d2f-8b63-ef623c9a778b	start-attP_Bxb1*-end
628bfe43-99c5-492e-8788-a52d387bc2a7	start-attP_Bxb1*-end
...
```

and a report of the run is printed to stdout. Its read accounting looks like this — see [Run report](#run-report) for the rest of it:

```
Summary
  reads read:                  100
  classified:                   62  (62.0%)
  not classified:               38  (38.0%)
    read too short              10
    read too long                8
    'start' segment not found    2
    'end' segment not found     17
    neither segment found        1
```

This should run very quickly and take less than a minute to run on a standard desktop machine.

Some further examples:

```sh
# Classify reads exactly as they arrive, with no anchoring or trimming.
segment --segments data/refs.fasta --sequences data/reads.fastq --classifications classifications.txt

# Also write a CSV counting how many reads produced each classification.
segment --segments data/refs.fasta --sequences data/reads.fastq -d \
        --classifications results.txt --counts counts.csv

# Take each segment's threshold from its name, as in ">SEG1[1.7]".
segment --segments scored_refs.fasta --sequences data/reads.fastq -d \
        --per-segment-scores --classifications classifications.txt

# Treat the construct as circular, using 8 threads and a stricter threshold.
segment --segments data/refs.fasta --sequences data/reads.fastq \
        --start-end-segs --circular --min-norm-score 1.7 --threads 8 \
        --classifications classifications.txt
```

## Input files

### Segments (`--segments`)

A FASTA file with one record per segment. The record ID (the first word of the header line) is the name used in the output, as in the `data/refs.fasta` file provided:

```
>start
CTCGGATACCCTTACTCTGTTGAAAACGAATAGATAGGTT
>end
ATTATTGACCACTTCCGAGTAGAATCGTGCTTCAGTAAGA

>attB_Tp901
ATGCCAACACAATTAACATCTCAATCAAGGTAAATGCTTTTTGCTTTTTTTGC
>attP_Tp901
GCGAGTTTTTATTTCGTTTATTTCAATTAAGGTAACTAAAAAACTCCTTT
...
```

Sequences are upper-cased on load, so a lower-case or soft-masked FASTA works fine. Segment names must be unique, and every record must have a sequence.

The names `start` and `end` are only special when `-d/--start-end-segs` is given. Without it they are ordinary segments, searched for and reported like any other.

#### Ambiguous bases

Segments may use the full set of IUPAC nucleotide codes, so a single segment can stand for several sequences at once — useful for a degenerate primer, a barcode with a wobble position, or a site that varies between variants:

| Code | Matches | Code | Matches |
|---|---|---|---|
| `A` `C` `G` `T` | themselves | `B` | C, G or T |
| `R` | A or G | `D` | A, G or T |
| `Y` | C or T | `H` | A, C or T |
| `S` | C or G | `V` | A, C or G |
| `W` | A or T | `N` | any base |
| `K` | G or T | | |
| `M` | A or C | | |

`U` is accepted and read as `T`, so an RNA spelling finds the same reads as the DNA one.

A code matches when the bases it permits overlap those the read permits, and a match scores the same `+2` whether the base was spelled out or reached through an ambiguity code. So `R` in a segment matches an `A` or a `G` in the read at full score, and costs an ordinary mismatch against a `C` or a `T`. The same rule applies in the other direction: an `N` a basecaller left in a read matches whatever the segment asks for at that position instead of costing a mismatch. Ambiguity codes work in the `start` and `end` anchors too, and are complemented correctly when segments are searched on the reverse strand (`R` becomes `Y`).

Because ambiguity is free under this scoring, a heavily degenerate segment can clear the score threshold on sequence that does not contain it — see [Choosing a score threshold](#choosing-a-score-threshold). A segment whose bases are *all* ambiguous would match at every position of every read, so it is refused outright, as is any segment containing a character that is not a nucleotide code at all.

#### Per-segment thresholds

One threshold rarely suits every segment: a long, distinctive one can afford a loose score, while a short or degenerate one needs a strict one to stay meaningful. With `--per-segment-scores`, a segment can name its own in square brackets after the record ID:

```
>start
CTCGGATACCCTTACTCTGTTGAAAACGAATAGATAGGTT
>attB_Tp901[1.7]
ATGCCAACACAATTAACATCTCAATCAAGGTAAATGCTTTTTGCTTTTTTTGC
>attP_Bxb1[1.2]
GGTTTGTCTGGTCAACCACCGCGGTCTCAGTGGTGTACGGTACAAACC
>attB_Int5
GAGCGCCGGATCAGGGAGTGGACGGCCTGGGAGCGCTACACGCTGTGGCTGCGGTCGGTGC
```

Here `attB_Tp901` is only reported at `1.7` or better and `attP_Bxb1` at `1.2` or better, while `start` and `attB_Int5` say nothing and so use `--min-norm-score`. The value is held to the same `0.0`–`2.0` range as `--min-norm-score`. Anchors are segments like any other, so `>start[1.7]` works and applies to that anchor alone.

**A bracketed value is always taken off the name**, whether or not the flag is given, so a segment is reported under the same name however the run is configured — the output reads `attB_Tp901`, never `attB_Tp901[1.7]`. The flag decides only whether the value inside is *used*. Without it the value is ignored, the segment is held to `--min-norm-score` like any other, and the run says so:

```
warning: segment 'attB_Tp901[1.7]' in 'refs.fasta' names a minimum normalised score, but
--per-segment-scores was not given, so it is ignored and the segment is held to
--min-norm-score like any other. The brackets are dropped from its name either way.
```

Because the brackets come off before names are compared, `SEG1[1.5]` and `SEG1[1.7]` are one segment asking for two thresholds, and are rejected as a duplicate.

#### When a score cannot be read

With the flag on, the value has to be a usable score. Brackets that cannot be read as one stop the run rather than being folded into the name, since a segment silently renamed is a segment silently never found. Every message names the record it came from — the ID as written, the file, and the record number — and says what specifically was wrong:

```
error: Segment 'attB_Tp901[high]' in 'refs.fasta' (record 3): 'high' is not a number.
With --per-segment-scores the square brackets hold that segment's minimum normalised
score, between 0.0 and 2.0, as in 'SEG1[1.5]'.
```

| Written as | Reported as |
|---|---|
| `SEG1[high]`, `SEG1[1.5.2]` | `'high' is not a number` |
| `SEG1[]` | `the square brackets are empty` |
| `SEG1[NaN]`, `SEG1[inf]` | `'NaN' is not a usable score` |
| `SEG1[2.5]` | `outside the range 0.0 to 2.0 … no alignment can score above 2.0, so the segment could never be found` |
| `SEG1[-0.5]` | `outside the range 0.0 to 2.0 … every alignment already scores at least 0.0, so the segment would be found everywhere` |
| `[1.5]` | `there is no name in front of the brackets` |
| `SEG1[1.5` | `the square bracket is never closed` |
| `SEG1[1.5]x` | `the score has to come last, but 'x' follows the closing bracket` |

Out-of-range scores say what would have happened rather than repeating the bounds back, and the value is quoted as it was written, so `[2.50]` is reported as `'2.50'` and not reformatted.

Only a *closed* bracket group at the end of the ID counts, so a name ending in `]` with no `[` before it — or one that simply contains brackets, like `SEG[1][1.5]` — keeps everything but the final pair. Without the flag nothing inside the brackets is interpreted at all, so `gene[human]` loads as the segment `gene` rather than being rejected.

### Reads (`--sequences`)

One of FASTA, FASTQ, either of those gzipped, or unaligned BAM (uBAM). The format is detected by inspecting the file contents — `>` for FASTA, `@` for FASTQ — so the extension does not matter and there is no flag to set. Because BAM is itself gzip-framed, detection looks past the shared gzip header at the decompressed content to tell a gzipped FASTA or FASTQ from a BAM.

FASTA records may wrap their sequence across as many lines as they like; each record is rejoined into one read. Qualities are never used for classification, so a FASTA run and a FASTQ run over the same sequences produce identical results.

Only *unaligned* BAM is accepted. An aligned BAM repeats a read once per alignment (secondary and supplementary records) and hard-clips supplementary sequence, which would silently inflate counts and truncate reads, so it is rejected with an error. Convert first:

```sh
samtools fastq aligned.bam > reads.fastq
```

Read sequences are upper-cased on load, and a uBAM record flagged as reverse complemented has its sequence flipped back into the orientation it was sequenced in.

## Output format

A tab-separated file with one line per accepted read and no header:

```
<read name><TAB><segment string>
```

[`--detailed-output`](#detailed-output) adds two more columns to each line.

The segment string lists the segments found, in the order they occur along the read, joined by `-`:

| String | Meaning |
|---|---|
| `attP_Bxb1-attB_Tp901` | Both segments found, in that order, on the forward strand |
| `attB_Tp901*` | Found in reverse-complement orientation (trailing `*`) |
| `attP_Bxb1-attP_Bxb1` | The same segment found twice |
| *(empty)* | No segment scored above the threshold |

With `-d/--start-end-segs`, the string is additionally bounded by the anchors, which are reported as the bounds rather than searched for within the read:

```
start-attP_Bxb1-attB_Tp901*-end
start-end                          ← anchors found, nothing in between
```

Note that reads which are filtered out do not appear in the output at all. A read is dropped if it falls outside the length bounds, or — with `-d` — if either anchor is missing, too degraded to score above the threshold, or in an inconsistent order. The [run report](#run-report) accounts for every one of them.

Results are deterministic: the same input, segments and settings always produce byte-identical output, including the order of lines, regardless of `--threads`.

While running, a progress bar on stderr shows how much of the file has been consumed, with a running read count and an ETA. It is measured in bytes of the input file, so it works the same for gzipped and BAM input, where the number of reads is not known until the file has been read.

### Detailed output

`--detailed-output` adds two columns to the same file:

```
<read name><TAB><segment string><TAB><located segments><TAB><extracted sequence>
```

```
ad6f445e   start-attP_Bxb1*-end   start[1:40],attP_Bxb1*[929:978],end[979:1018]   CTCGGATACCCTTA…
```

**Located segments** repeats the segment string with the span each hit covers, and does so name for name: every name in the segment string has an entry here, in the same order. Entries are separated by commas, and the two positions within an entry by a colon, so the column can be split on `,` without the two separators being confused for one another.

Positions are **1-based and inclusive of both ends**, counted along the extracted sequence in the next column. So `attP_Bxb1*[929:978]` is 50 bases and `sequence[929..=978]` cuts it out. Reverse-strand hits keep the same trailing `*` they have in the segment string. A read with no segments above the threshold gets an empty column rather than a missing one, so every row has the same number of columns.

With `-d/--start-end-segs` the anchors are reported here too. They are always the two ends of the extracted sequence — `start` opens it and `end` closes it — but the *spans* are worth having: they are what the anchors actually aligned across, so an anchor that matched short or long shows up here rather than having to be inferred. In the example above `start[1:40]` covers all 40 bases of the anchor, while a read whose anchor aligned with a deletion would report 39.

**Extracted sequence** is the read exactly as it was classified — the whole read when no anchors were given, and the trimmed, reoriented span from `start` to `end` (anchors included) when they were. This is what makes the positions well defined: with `-d` a read sequenced backwards is reverse complemented before classification, so both strands of the same molecule produce identical positions, counted along the sequence in the column rather than along the read as it arrived.

Because the extracted sequence is carried alongside each read until its chunk is written, this roughly doubles what a chunk holds in memory. Chunks are bounded, so the ceiling does not grow with the size of the input.

### Classification counts

`--counts <FILE>` writes a second file alongside the per-read results: a CSV counting how many reads produced each distinct segment string. It is an addition, never a substitute — the results file is written exactly as it would have been without the flag.

```sh
segment --segments refs.fasta --sequences reads.fastq -d \
        --classifications results.txt --counts counts.csv
```

```csv
segments,count
start-attP_Bxb1*-end,53
start-attB_Tp901-attB_Int5-end,7
start-attB_Tp901-attB_Int5-attB_BxB1*-attP_Int5*-attP_Tp901-attP_Bxb1*-end,5
start-end,4
start-attB_Int5-attB_BxB1*-attP_Int5*-attP_Tp901-attP_Bxb1*-end,1
```

Rows are ordered by count, most frequent first, with equal counts ordered alphabetically so that two runs over the same reads produce byte-identical files and a diff between two conditions stays readable.

The counts cover the reads that were *classified*, and so add up to the number of lines in the results file rather than to the number of reads in the input. Reads dropped for their length or for a missing anchor never reach a classification and are accounted for in the [run report](#run-report) instead. A read that was classified but carried no recognisable segments is a real result and gets a row with an empty first field.

A run that classified nothing still writes the header row, so the file is always valid CSV. Fields are quoted where a segment name contains a comma, a quote or a newline.

## How it works

Each read goes through up to four stages. Reads are streamed through in bounded chunks rather than loaded all at once, so a file larger than memory is fine. Chunks are classified and written in order, which is what keeps the output identical however large the input.

### 1. Length filter

Reads shorter than `--min-seq-len` or longer than `--max-seq-len` are dropped before any alignment work is done. Each bound is `0` by default, meaning no bound on that side, so out of the box nothing is dropped for its length. The two are independent: `--min-seq-len 500` on its own drops short reads and leaves the upper end unbounded.

### 2. Anchoring and orientation (only with `-d`)

A long read usually contains more than your construct: adapters, barcodes and other flanking sequence, and it may have been sequenced off either strand.

If you pass `-d/--start-end-segs`, the segments named `start` and `end` are used as anchors. Both are aligned against the read in both orientations, and the read is:

1. **rejected** if either anchor fails to score above the threshold;
2. **reverse complemented** if it was sequenced off the opposite strand, so that downstream results are strand-independent;
3. **trimmed** to the region running from the `start` anchor to the `end` anchor.

Everything outside the anchors is discarded and cannot contribute segments. A read whose `start` does not precede its `end` is inconsistent for a linear construct and is rejected — unless `-c/--circular` is set, in which case it is treated as crossing the origin and rotated back into order (see [below](#circular-constructs)).

Without `-d`, this stage is skipped entirely: the read is classified exactly as it arrives, on whichever strand it arrived.

### 3. Segment alignment

Every segment is aligned against the read in both orientations using a semi-global alignment — the whole segment must align, but it may sit anywhere in the read, and gaps at either end of the read are free. Scoring is `+2` for a match, `-1` for a mismatch, and `-2` for a gap (insertion or deletion). Two bases match when the sets of nucleotides they stand for overlap, which for plain `ACGT` is ordinary equality and for [ambiguous bases](#ambiguous-bases) is rather more forgiving.

The **normalised score** is the raw score divided by the segment length, so a perfect match scores `2.0` regardless of how long the segment is. A segment is only considered found where its normalised score reaches `--min-norm-score`.

Crucially, *every* position where a segment scores above the threshold is recorded, not just the best one — so a segment occurring three times in a read is found three times.

### 4. Overlap resolution

Because all segments are searched independently, their hits can overlap. These are resolved greedily, best-scoring first:

- Hits overlapping by **3 bp or less** are both kept. Segments in real constructs frequently abut or share a few bases, and this tolerance stops that being penalised.
- Hits overlapping by **more than 3 bp** cannot both be real, so the lower-scoring one is discarded. Because resolution runs in descending score order, the better alignment always wins regardless of where it sits in the read or where it appeared in the FASTA.
- A hit entirely **containing** a better-scoring hit is discarded in favour of the contained one.

Ties are broken by position and then by segment name, so two equally good alignments always resolve the same way. Surviving hits are then sorted by position and joined to produce the segment string.

### Circular constructs

A read from a circular template (a plasmid, say) can begin anywhere, so it may run `end`-…-`start` rather than `start`-…-`end`. With `-c/--circular` such a read is rotated about its origin back into `start`-…-`end` order before classification, on either strand. Without `-c` it is rejected as inconsistent. `-c` only affects reads that actually wrap; reads that do not cross the origin follow the ordinary path.

## Choosing a score threshold

`--min-norm-score` is the single most important parameter for result quality.

A normalised score of `s` on an ungapped segment needs a fraction `(s + 1) / 3` of the bases to match — so `1.3` needs about 77% identity and `1.5` about 83%. (Rounding up to a whole number of bases makes short segments slightly stricter: a 20 bp segment at `1.5` needs 17/20, or 85%.)

The looser the threshold, the more spurious hits appear, and the effect grows with read length because there are simply more places for a chance match to occur. Searching six 20 bp segments against 50 kb of *random* sequence containing none of them:

| `--min-norm-score` | Spurious segments reported |
|---|---|
| `1.7` | 0 |
| `1.5` (default) | 0 |
| `1.3` | 6 |
| `1.1` | 144 |
| `0.9` | 1038 |

The default of `1.5` reports nothing on the random control above while still tolerating about one mismatch in six, which suits noisy long reads. Lower it only if you have reason to — very noisy data, or segments long enough that a few percent of chance identity cannot accumulate into a false hit — and whenever you do, sanity-check by running against sequence you know does not contain your segments.

### Setting a threshold per segment

`--min-norm-score` applies to every segment at once, which is a blunt instrument when segments differ in length and specificity. [`--per-segment-scores`](#per-segment-thresholds) lets each one name its own in its FASTA header, so a short or degenerate segment can be held to a strict score without making every other segment strict too. Everything below applies per segment rather than globally once that flag is on.

### Thresholds and ambiguous bases

[Ambiguity codes](#ambiguous-bases) push the same effect the other way. A spelled-out base matches one random base in four, but an `R` matches one in two and an `N` matches every one, so a degenerate segment finds chance matches far more readily than its length suggests. The identity a threshold demands is fixed — `(s + 1) / 3`, so 83% at the default `1.5` — while the identity a segment achieves by luck alone is the average of `bases permitted / 4` across its positions.

When the second figure reaches the first, the segment clears the threshold on random sequence and will be reported almost everywhere. `segment` checks this on startup and names any segment that qualifies:

```
warning: segment 'wobble_probe' is ambiguous enough to match 85% of random bases, which
is at or above the 83% identity its minimum normalised score asks for, so expect to find
it almost anywhere. Use a more specific sequence, or raise the score it is found at.
```

This is a warning rather than an error, since a mostly degenerate probe is a legitimate thing to search for. Each segment is judged against its own threshold, so giving just that segment a stricter one with [`--per-segment-scores`](#per-segment-thresholds) both silences the warning and fixes what it is warning about, without making every other segment stricter.

## Run report

At the end of every run a report is written to **stdout**, describing what went in, how the run was configured, how the reads were accounted for, and how fast it went:

```
Run
  started                      2026-09-01 18:24:13 UTC
  finished                     2026-09-01 18:24:13 UTC
  elapsed                      0.31 s

Input
  segments file                data/refs.fasta
  segments loaded                8
  sequences file               data/reads.fastq
  sequences format             FASTQ
  sequences size               270.4 KiB

Options
  --min-norm-score             1.5
  --per-segment-scores         off
  --start-end-segs             on
  --circular                   off
  --min-seq-len                0
  --max-seq-len                0
  --threads                    1
  --detailed-output            off
  --classifications            classifications.txt
  --counts                     (not written)

Summary
  reads read:                  100
  classified:                   70  (70.0%)
  not classified:               30  (30.0%)
    'start' segment not found    2
    'end' segment not found     24
    neither segment found        4

Throughput
  bases read                   126.0 kbase
  reads per second             321
  bases per second             404.0 kbase
```

**Run** timestamps the run. Times are UTC so they sort and compare without ambiguity about the machine's timezone. `elapsed` is measured on a monotonic clock, so an adjustment to the system clock mid-run cannot distort it.

**Input** names both files, how many segments were loaded from the first, and the detected format and size of the second — so a report says which files produced it, not just what came out.

**Options** lists *every* option at the value it actually ran with, including the ones left at their defaults. Together with the input section this makes the report a record of what produced the results beside it, which is what makes a run reproducible from its own output.

**Summary** accounts for every read. Classified and not-classified always add back up to the number of reads, and each rejected read is counted under exactly one reason. Categories that did not occur are left out, so a run with no anchoring shows only the length rejections.

**Throughput** is the rate over the whole run, loading included. A run too short for the clock to measure quotes no rate rather than an infinity.

Every count shares one column whatever its nesting depth, and the columns are measured from the rows actually being printed, so a run of a hundred reads and a run of ten million both come out square:

```
Summary
  reads read:      100
  classified:      100  (100.0%)
  not classified:    0    (0.0%)
```

| Reason | Meaning |
|---|---|
| `read too short` / `read too long` | Outside `--min-seq-len` / `--max-seq-len`. Counted before any alignment. |
| `'start' segment not found` | The `start` anchor scored below the threshold on both strands. Usually a truncated read or a threshold set too high. |
| `'end' segment not found` | As above, for the `end` anchor. |
| `neither segment found` | Neither anchor was found. Often an off-target read, or the wrong segments file. |
| `'start'/'end' on opposite strands` | Both anchors were found, but facing opposite ways, so the read cannot be oriented. Typically a chimeric read. |
| `'start'/'end' out of order` | Both anchors found on the same strand, but `end` precedes `start`. Add `--circular` if the template is circular and these reads cross the origin. |
| `records skipped while reading the input` | Malformed, unnamed or sequence-less records skipped during loading. Each is also warned about individually. |

Anchor-related reasons only arise with `-d/--start-end-segs`.

The report goes to stdout while warnings go to stderr, and the results go to their own file, so the three can be captured separately. Nothing else is written to stdout, so `> run.txt` keeps the report and nothing but the report:

```sh
segment --segments refs.fasta --sequences reads.fastq -d \
        --classifications results.txt > run.txt 2> warnings.log
```

## Errors and warnings

`segment` distinguishes between problems that invalidate the whole run and problems affecting a single record.

**Errors** stop the run, print a message to stderr, and exit with status `1`. These are setup mistakes — a missing or unreadable file, a `--classifications` or `--counts` path that cannot be created, a segments file with duplicate names, an empty record, a character that is not a nucleotide code, a wholly ambiguous segment, an unreadable per-segment score, an aligned BAM, `-d` without `start`/`end` segments, or an out-of-range argument:

```
error: Segment 'attB_Tp901' is defined more than once in 'parts.fasta'. Segment names
must be unique, otherwise it is ambiguous which sequence should be searched for.
```

**Warnings** go to stderr and processing continues. On startup, a segment [degenerate enough to match random sequence](#thresholds-and-ambiguous-bases) is named, as is one that [asks for a threshold](#per-segment-thresholds) without `--per-segment-scores` to apply it. An individual malformed or unusable record is skipped, named where possible, and a count of them is printed once the file is read:

```
warning: skipping malformed read 'ad6f445e' in 'reads.fastq': Incomplete record. ...
warning: skipping read '124e2e64' in 'reads.fastq': it has no sequence
warning: skipped 2 unusable record(s) in 'reads.fastq'; 41998 read(s) loaded
```

One truncated record at the end of a large FASTQ therefore costs you that record, not the run. A BAM whose *compressed stream* is corrupt is a fatal error rather than a warning, because there is no safe way to resynchronise and continue.

## Development

Run the test suite:

```sh
cargo test
```

The suite is self-contained — every fixture (FASTA, FASTQ, gzipped FASTQ and BAM) is generated into a temporary directory at test time, so no test data is kept in the repository. It covers input detection and loading, alignment scoring, overlap resolution, anchoring and reorientation, circular constructs, input validation, and end-to-end runs of the compiled binary.

Two earlier implementations are kept purely as reference oracles, so the rewrites that replaced them stay honest:

- The original full-matrix segment aligner, in `src/tests.rs`. It is compared against the current linear-space one exhaustively over every short sequence pair from a two-letter alphabet, over randomised inputs and at realistic scale, checking every field of every alignment agrees.
- The library aligner previously used for the `start`/`end` anchors. The replacement reproduces its affine gap model, and scores — which decide whether a read is accepted — are checked for exact equality. On deliberately degenerate sequence the two can report different, equally-scoring alignments, since which of several co-optimal alignments to report is an arbitrary choice; on realistic reads they agree completely.

Layout:

| Path | Contents |
|---|---|
| `src/main.rs` | The tool: I/O, the aligner, classification, and the CLI. |
| `src/tests.rs` | Unit tests, with access to the internals. |
| `tests/cli.rs` | End-to-end tests driving the compiled binary. |

## Licence

See [LICENSE](LICENSE).
