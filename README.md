> Warning: This is a pre-release version, we are actively developing this repository. Issues, bugs and features will happen, rise and change.

# COVIDmutants

We foster the openness, integrity, and reproducibility of scientific research.

## How to use this repository?

This repository host both the scripts and tools developed by this study. Feel free to adapt the scripts and tools, but remember to cite their authors!

To look at our scripts and raw results, **browse** through this repository. If you want to reproduce our results you will need to **clone** this repository, build the docker, and the run all the scripts. If you want to use our data for our own research, **fork** this repository and **cite** the authors.


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
  -o "$(pwd)/test"
  -fasta "$(pwd)/resources/gisaid_hcov-19_2021_02_12_21.fasta.gz" \
```

#### Start the analysis

With a `covid_genomes.fasta` as your multi-sequence fasta files (gzip compressed file are accepted).

```sh
docker run -i -t --rm -v "$(pwd)":"$(pwd)" -u $(id -u):$(id -g) covidmutants \
  -threads 15 \
  -forward TATgCCATTAgTgCAAAgAATAgAgCTCgCAC \
  -reverse GTAATTGGAAcAAGcAAATTcTATGGTGGTTG \
  -amplicon TATGCCATTAGTGCAAAGAATAGAGCTCGCACCGTAGCTGGTGTCTCTATCTGTAGTACTATGACCAATAGACAGTTTCATCAAAAATTATTGAAATCAATAGCCGCCACTAGAGGAGCTACTGTAGTAATTGGAACAAGCAAATTCTATGGTGGTTG \
  -o "$(pwd)/output"
  -fasta "$(pwd)/covid_genomes.fasta(.gz)" \
```

## Usage

```plaintext
usage: docker run -i -t --rm -v $(pwd):$(pwd) covidmutants
                [-h] -fasta GENOMESFILE -forward FORWARD -reverse REVERSE
                [-probe PROBE] [-amplicon AMPLICON] [-output OUTFILE]
                [-threads THREADS] [-all] [-json]

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
  -output OUTFILE, -out OUTFILE, -o OUTFILE
                        Specify output file
  -threads THREADS, -t THREADS
                        Maximum threads
  -all, -a              Test all sequences (VERY slow)
  -json                 Export results matrix as JSON file (Experimental)
```

## Methodology

...

## Results

...

```plaintext
GENOME_ID       PASS/FAIL   FORWARD_DIFF                          REVERSE_DIFF                        AMPLICON_DIFF
ZMB-V65185      PASS        ................................     ................................
NB-RIVM-21079   PASS        ................................     ................................
CAMC-11B5E57    PASS        ................................     ................................
N00641/2020     FAIL        ........A......CCt.c.c....t.....     ................................     .......................................................
N00663/2020     PASS        ................................     ................................
IDF-IPP01211    PASS        ................................     ................................
```

## Issues

If you have any problems with or questions about the scripts, please contact us through a [GitHub issue](https://github.com/aicuramedical/covidmutants/issues).
Any issue related to the scientific results themselves must be done directly with the authors.


## Contributing

You are invited to contribute new features, fixes, or updates, large or small; we are always thrilled to receive pull requests, and do our best to process them as fast as we can.


## License and distribution

The content of this project itself including the raw data and work are licensed under the [Creative Commons Attribution-ShareAlike 4.0 International License](http://creativecommons.org/licenses/by-sa/4.0/), and the source code presented is licensed under the [GPLv3 license](http://www.gnu.org/licenses/gpl-3.0.html).







