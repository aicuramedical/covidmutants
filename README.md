> Warning: This is a pre-release version, we are actively developing this repository. Issues, bugs and features will happen, rise and change.

# COVIDmutants

We foster the openness, integrity, and reproducibility of scientific research.


## How to use this repository?

This repository hosts the scripts and tools created for research.
Please feel free to modify the scripts, please remember to credit the creators!


## Prepare a docker

All required files and tools run in a self-contained [docker](https://www.docker.com/) image.

#### Clone the repository

```sh
git clone https://github.com/aicuramedical/covidmutants.git
cd covidmutants
```

#### Create a docker

```sh
docker build --rm=true -t covidmutants .
```

#### Run a test

```sh
docker run -i -t --rm -v "$(pwd)":"$(pwd)" covidmutants \
  -threads 15 \
  -forward TATgCCATTAgTgCAAAgAATAgAgCTCgCAC \
  -reverse GTAATTGGAAcAAGcAAATTcTATGGTGGTTG \
  -amplicon TATGCCATTAGTGCAAAGAATAGAGCTCGCACCGTAGCTGGTGTCTCTATCTGTAGTACTATGACCAATAGACAGTTTCATCAAAAATTATTGAAATCAATAGCCGCCACTAGAGGAGCTACTGTAGTAATTGGAACAAGCAAATTCTATGGTGGTTG \
  -o "$(pwd)"/test \
  -v \
  -fasta "$(pwd)"/test/test.fa.gz
```

#### Start the analysis

With a `covid_genomes.fasta` as your multi-sequence fasta files (gzip compressed file are accepted).

```sh
docker run -i -t --rm -v "$(pwd)":"$(pwd)" -u $(id -u):$(id -g) covidmutants \
  -threads 15 \
  -forward TATgCCATTAgTgCAAAgAATAgAgCTCgCAC \
  -reverse GTAATTGGAAcAAGcAAATTcTATGGTGGTTG \
  -amplicon TATGCCATTAGTGCAAAGAATAGAGCTCGCACCGTAGCTGGTGTCTCTATCTGTAGTACTATGACCAATAGACAGTTTCATCAAAAATTATTGAAATCAATAGCCGCCACTAGAGGAGCTACTGTAGTAATTGGAACAAGCAAATTCTATGGTGGTTG \
  -o "$(pwd)"/output \
  -fasta "$(pwd)"/covid_genomes.fasta(.gz)
```


## Usage

```plaintext
usage: docker run -i -t --rm -v $(pwd):$(pwd) covidmutants
                [-h] -fasta GENOMESFILE -forward FORWARD -reverse REVERSE
                [-probe PROBE] [-amplicon AMPLICON] [-max MAX_DIFF]
                [-output OUTFILE] [-threads THREADS] [-all] [-json]


optional arguments:
  -h, --help            show this help message and exit
  -fasta GENOMESFILE, -i GENOMESFILE, -in GENOMESFILE
                        Specify covid fasta genomes file
  -forward FORWARD, -1 FORWARD
                        Specify forward primer
  -reverse REVERSE, -2 REVERSE
                        Specify reverse primer
  -probe PROBE          Specify probe sequence
  -amplicon AMPLICON    Specify amplicon sequence
  -max MAX_DIFF, -m MAX_DIFF
                        Maximum sequence variations before querying sequencing accuracy
  -output OUTFILE, -out OUTFILE, -o OUTFILE
                        Specify output file
  -threads THREADS, -t THREADS
                        Maximum threads
  -all, -a              Test all sequences (VERY slow)
  -json                 Export results matrix as JSON file (Experimental)
```


## Methodology

[to be improved]

1. Validate and clean up the sequences (check name and duplication)
2. Use KAT to quickly recover the exact match (k-mer = primer length).
3. Use Smith-Waterman local alignment to find sequence variations.
4. Report


## Results

#### \<output\>.tsv file

For each sequence submitted, the name and variation with the forward and reverse primers are reported. If the amplicon and/or probe sequence were provided, the sequence variations are reported only if the primers failed to pass the filter (or the parameter `-all` was specified.)


```plaintext
GENOME_ID       PASS/FAIL   FORWARD_DIFF                          REVERSE_DIFF                        AMPLICON_DIFF
ZMB-V65185      PASS        ................................     ................................
NB-RIVM-21079   PASS        ................................     ................................
CAMC-11B5E57    PASS        ................................     ................................
N00641/2020     FAIL        ........A......CCt.c.c....t.....     ................................     .......................................................
N00663/2020     PASS        ................................     ................................
IDF-IPP01211    PASS        ................................     ................................
```

#### \<output\>.fail.fa file

FAIL indicates that the primer sequences did not match the primers submitted. An extra file \<output\>.fail.fa is created with those sequences for manual validation.

#### \<output\>.issue.fa file

If the sequences (primer or amplicon) are found but have too many sequence variants (e.g. more than `-max`), the accuracy of their sequences is questioned.
The sequences are reported separately in the \<output\>.issue.fa file.


## Issues

If you have any problems with or questions about the scripts, please contact us through a [GitHub issue](https://github.com/aicuramedical/covidmutants/issues).
Any issue related to the scientific results themselves must be done directly with the authors.


## Contributing

You are invited to contribute new features, fixes, or updates, large or small; we are always thrilled to receive pull requests, and do our best to process them as fast as we can.


## License and distribution

The content of this project itself including the raw data and work are licensed under the [Creative Commons Attribution-ShareAlike 4.0 International License](http://creativecommons.org/licenses/by-sa/4.0/), and the source code presented is licensed under the [GPLv3 license](http://www.gnu.org/licenses/gpl-3.0.html).







