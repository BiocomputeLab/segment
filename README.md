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
segment --segments REF_FASTA --sequences READS [OPTIONS] --output OUTPUT_FILE
```

where `REF_FASTA` is a file containing the reference segments in FASTA format, `READS` is the sequencing data (FASTQ, gzipped FASTQ, or unaligned BAM), and `OUTPUT_FILE` is the filename to write the results to.

By default every read is classified exactly as it arrives. Passing `-d/--start-end-segs` instead uses the segments named `start` and `end` to anchor each read: reads are oriented, trimmed to the region between the two anchors, and dropped if either anchor is missing. This is usually what you want when the reads contain adapters or other flanking sequence.

### Required arguments

| Argument | Description |
|---|---|
| `--segments <FILE>` | FASTA file of segment sequences to search for, one record per segment. Record IDs become the names used in the output. Names must be unique and no record may be empty. |
| `--sequences <FILE>` | Reads to classify: FASTQ, gzipped FASTQ, or unaligned BAM. The format is detected from the file contents, so the extension is irrelevant. |

### Optional arguments

| Argument | Default | Description |
|---|---|---|
| `-o`, `--output <FILE>` | `output.txt` | Where to write the tab-separated results. Overwritten if it exists. |
| `-d`, `--start-end-segs` | off | Use the segments named `start` and `end` to anchor, trim and orient each read; they then bound the segment string instead of being reported within it. Reads missing either anchor are dropped. Errors if the segments file has no `start` or `end` record. |
| `-c`, `--circular` | off | Treat reads as coming from a circular template, so a read that crosses the origin (`end` before `start`) is rotated back into order rather than rejected. Only meaningful together with `-d`. |
| `-s`, `--min-norm-score <FLOAT>` | `1.5` | Minimum normalised alignment score for a segment to count as found, from `0.0` to `2.0`, where `2.0` is a perfect match. See [Choosing a score threshold](#choosing-a-score-threshold). |
| `--min-seq-len <INT>` | `100` | Reads shorter than this are dropped before alignment. Set this and `--max-seq-len` both to `0` to disable length filtering. |
| `--max-seq-len <INT>` | `100000` | Reads longer than this are dropped before alignment. Must not be less than `--min-seq-len`. |
| `-t`, `--threads <INT>` | `1` | Reads are classified in parallel across this many threads. |
| `-h`, `--help` | | Print help. |
| `-V`, `--version` | | Print version. |

Peak memory does not grow with the size of the input. Reads are streamed through in bounded chunks, and all alignment runs in space linear in the read length and independent of segment length, so long segments and long anchors cost no more than short ones. The practical ceiling is a single very long read, not the file.

## Example

To run `segment` on the small example data set provided in the `data` directory, we would therefore use the command (running from the root of the code repository):

```sh
segment --segments data/refs.fasta --sequences data/reads.fastq \
        --start-end-segs --min-seq-len 500 --max-seq-len 2000 \
        --output output.txt
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

and a summary of the run is printed to stderr:

```
Summary
  reads read:            100
  classified:             62  (62.0%)
  not classified:         38  (38.0%)
    read too short                           10
    read too long                             8
    'start' segment not found                 2
    'end' segment not found                  17
    neither segment found                     1
```

This should run very quickly and take less than a minute to run on a standard desktop machine.

Some further examples:

```sh
# Classify reads exactly as they arrive, with no anchoring or trimming.
segment --segments data/refs.fasta --sequences data/reads.fastq --output output.txt

# Treat the construct as circular, using 8 threads and a stricter threshold.
segment --segments data/refs.fasta --sequences data/reads.fastq \
        --start-end-segs --circular --min-norm-score 1.7 --threads 8 \
        --output output.txt
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

### Reads (`--sequences`)

One of FASTQ, gzipped FASTQ, or unaligned BAM (uBAM). The format is detected by inspecting the file contents, so the extension does not matter and there is no flag to set. Because BAM is itself gzip-framed, detection looks past the shared gzip header at the decompressed content to tell a gzipped FASTQ from a BAM.

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

Note that reads which are filtered out do not appear in the output at all. A read is dropped if it falls outside the length bounds, or — with `-d` — if either anchor is missing, too degraded to score above the threshold, or in an inconsistent order. The [run summary](#run-summary) accounts for every one of them.

Results are deterministic: the same input, segments and settings always produce byte-identical output, including the order of lines, regardless of `--threads`.

While running, a progress bar on stderr shows how much of the file has been consumed, with a running read count and an ETA. It is measured in bytes of the input file, so it works the same for gzipped and BAM input, where the number of reads is not known until the file has been read.

## How it works

Each read goes through up to four stages. Reads are streamed through in bounded chunks rather than loaded all at once, so a file larger than memory is fine. Chunks are classified and written in order, which is what keeps the output identical however large the input.

### 1. Length filter

Reads shorter than `--min-seq-len` or longer than `--max-seq-len` are dropped before any alignment work is done. Setting both to `0` disables the check.

### 2. Anchoring and orientation (only with `-d`)

A long read usually contains more than your construct: adapters, barcodes and other flanking sequence, and it may have been sequenced off either strand.

If you pass `-d/--start-end-segs`, the segments named `start` and `end` are used as anchors. Both are aligned against the read in both orientations, and the read is:

1. **rejected** if either anchor fails to score above the threshold;
2. **reverse complemented** if it was sequenced off the opposite strand, so that downstream results are strand-independent;
3. **trimmed** to the region running from the `start` anchor to the `end` anchor.

Everything outside the anchors is discarded and cannot contribute segments. A read whose `start` does not precede its `end` is inconsistent for a linear construct and is rejected — unless `-c/--circular` is set, in which case it is treated as crossing the origin and rotated back into order (see [below](#circular-constructs)).

Without `-d`, this stage is skipped entirely: the read is classified exactly as it arrives, on whichever strand it arrived.

### 3. Segment alignment

Every segment is aligned against the read in both orientations using a semi-global alignment — the whole segment must align, but it may sit anywhere in the read, and gaps at either end of the read are free. Scoring is `+2` for a match, `-1` for a mismatch, and `-2` for a gap (insertion or deletion).

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

## Run summary

At the end of every run a summary is written to stderr, accounting for every read:

```
Summary
  reads read:             10
  classified:              4  (40.0%)
  not classified:          6  (60.0%)
    read too short                            1
    'start' segment not found                 1
    'end' segment not found                   1
    neither segment found                     1
    'start'/'end' on opposite strands         1
    'start'/'end' out of order                1
  records skipped while reading the input: 1
```

Classified and not-classified always add back up to the number of reads, and each rejected read is counted under exactly one reason. Categories that did not occur are left out, so a run with no anchoring shows only the length rejections.

| Reason | Meaning |
|---|---|
| `read too short` / `read too long` | Outside `--min-seq-len` / `--max-seq-len`. Counted before any alignment. |
| `'start' segment not found` | The `start` anchor scored below the threshold on both strands. Usually a truncated read or a threshold set too high. |
| `'end' segment not found` | As above, for the `end` anchor. |
| `neither segment found` | Neither anchor was found. Often an off-target read, or the wrong segments file. |
| `'start'/'end' on opposite strands` | Both anchors were found, but facing opposite ways, so the read cannot be oriented. Typically a chimeric read. |
| `'start'/'end' out of order` | Both anchors found on the same strand, but `end` precedes `start`. Add `--circular` if the template is circular and these reads cross the origin. |
| `records skipped while reading the input` | Malformed, unnamed or sequence-less records skipped during loading. Each is also warned about individually. |

Anchor-related reasons only arise with `-d/--start-end-segs`. Because the summary goes to stderr it never mixes with the results file, and `2>run.log` captures it together with any warnings.

## Errors and warnings

`segment` distinguishes between problems that invalidate the whole run and problems affecting a single record.

**Errors** stop the run, print a message to stderr, and exit with status `1`. These are setup mistakes — a missing or unreadable file, a segments file with duplicate names or an empty record, an aligned BAM, `-d` without `start`/`end` segments, or an out-of-range argument:

```
error: Segment 'attB_Tp901' is defined more than once in 'parts.fasta'. Segment names
must be unique, otherwise it is ambiguous which sequence should be searched for.
```

**Warnings** go to stderr and processing continues. An individual malformed or unusable record is skipped, named where possible, and a summary is printed once the file is read:

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
