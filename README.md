# segment

A flexible analysis tool for classifying the order and orientation of sequence segments in recombinase-based genetic circuitry.

## Installation

To use `segment` you will need to first compile the source code contained in this repository. This requires that the `rustc` compiler is available on your system. Further details about how to install `rustc` can be found at the following URL and we recommend the use of `rustup` to maintain your environment and ensure all other supporting tools are available.

https://www.rust-lang.org/tools/install

Once `rustc` is installed, compilation of `segment` is performed by running the following command from within the root of the code repository:

``sh
cargo build --release
``

This should take a couple minutes to complete and the generated `segment` executable can be found in the `target/release` directory. We recommend placing the `segment` executable in a location found in your `PATH` environment variable to making running the tool easier from the command line.

It should be possible to compile and run `segment` on any system for which `rustc` is available (e.g., Windows, MacOS and Linux) and the only dependancy is on the `bio` crate. Testing of the software has been performed exclusively on MacOS Sequoia (15.5).

## Usage

To run `segment` it has the following usage:

``sh
segment REF_FASTA DATA_FASTQ MIN_READ_LEN MAX_READ_LEN OUTPUT_FILE
``

where `REF_FASTA` is the location of a file containing the reference segments in a FASTA format (this must include a `start` and `end` segment reference as a bare minimum), `DATA_FASTQ` is the sequencing data in a FASTQ format, `MIN_READ_LEN` and `MAX_READ_LEN` are the minimum and maximum read length, respectively, to carry out the classification for, and `OUTPUT_FILE` is the filename to output the results to.

To run `segment` on the small example data set provided in the `data` directory, we would therefore use the command (running from the root of the code repository):

``sh
segment data/refs.fasta data/reads.fastq 500 2000 output.txt
``

This should generate a file classifying each read in the intput FASTQ data file that passes the read length filtering. The first few lines are:

``
ad6f445e-17e6-49ce-bf40-67da322006c8    start-attP_Bxb1*-end
124e2e64-ee6f-4349-972d-58e5a210c264    start-attP_Bxb1*-end
27bcf6c4-9226-4fa9-af0a-1ed4294d8768    start-end
715b0122-721e-4dc1-acff-ca723162392c    start-attP_Bxb1*-end
5f2bb513-496e-447b-a151-df14796f39bc    start-attP_Bxb1*-end
fe1b90af-4e41-4c33-ab87-eee17512462c    start-attP_Bxb1*-end
cb24e39a-8ee8-4988-a788-250ccd99349c    start-attB_Tp901-attB_Int5-attB_BxB1*-end
c689dd31-4e98-4154-ab05-a8b0abed9436    start-attP_Bxb1*-end
8c19bcd4-b228-4e0d-a257-e5116a0a4322    start-attP_Bxb1*-end
``
