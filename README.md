# COVIDmutants
COVIDmutants is a fast and scalable pipeline to check for the intact presence of SARS-CoV-2 Primers in an enormous collection of genome sequences and to report the findings.

## Associated publications

> **Efficient screening of long oligonucleotides against hundred thousands of SARS-CoV-2 genome sequences**.
> Weidmann M, Graf E, Lichterfeld D, El Wahed A and Bekaert M.
> _Front.Virol._ 2022

[![DOI](https://img.shields.io/badge/DOI-10.3389%2Ffviro.2022.835707-blue.svg)](https://doi.org/10.3389/fviro.2022.835707)

## How to use this repository?

We foster the openness, integrity, and reproducibility of scientific research. This repository hosts the scripts and tools created for research.
Please feel free to modify the scripts, please remember to credit the authors.

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
docker run -i -t --rm covidmutants \
  -threads 2 \
  -forward TATgCCATTAgTgCAAAgAATAgAgCTCgCAC \
  -reverse GTAATTGGAAcAAGcAAATTcTATGGTGGTTG \
  -amplicon TATGCCATTAGTGCAAAGAATAGAGCTCGCACCGTAGCTGGTGTCTCTATCTGTAGTACTATGACCAATAGACAGTTTCATCAAAAATTATTGAAATCAATAGCCGCCACTAGAGGAGCTACTGTAGTAATTGGAACAAGCAAATTCTATGGTGGTTG \
  -fasta test.fa.gz
```

## Dependencies

If you use the docker, we will not need to check dependencies, as the Dockerfile will contain all you need. But here are the details anyway:

* [KAT](https://github.com/TGAC/KAT) v2.4.2+
* [EMBOSS](http://emboss.open-bio.org/) v6.6.0+ suit, only [_water_](http://emboss.open-bio.org/rel/rel6/apps/water.html) is required for local alignments
* [Python](https://www.python.org/) v3.11+ with [biopython](https://biopython.org/)

## Methodology

We used a fast and scalable approach in order to screen the ever-growing number of genome sequences. To do so, we used K-mer of the exact size of the primer to mind through the genome using KAT mode "sect" ([Mapleson 2017](https://doi.org/10.1093/bioinformatics/btw663)). Kat calculates the K-mer coverage across every genome and reported the primer-coverage/integrity. Sequence that reported imperfect primer presence, were then properly aligned against the primers and the amplicon using Smith-Waterman ([Smith 1981](https://doi.org/10.1016/0022-2836%2881%2990087-5)) local alignment to find sequence variations. A custom python script handles the pipeline and reporting.

## Usage

```plaintext
usage: docker run -i -t --rm -v $(pwd):$(pwd) covidmutants
                [-h] -fasta GENOMESFILE -forward FORWARD -reverse REVERSE
                [-probe PROBE] [-amplicon AMPLICON] [-max MAX_DIFF]
                [-output OUTFILE] [-threads THREADS] [-all] [-json]


  -h, --help            show this help message and exit
  -fasta GENOMESFILE, -i GENOMESFILE, -in GENOMESFILE
                        Specify covid fasta genomes file
  -forward FORWARD, -1 FORWARD
                        Specify forward primer
  -reverse REVERSE, -2 REVERSE
                        Specify reverse primer
  -probe PROBE          Specify probe sequence(s) [can be used multiple times]
  -amplicon AMPLICON    Specify amplicon sequence
  -max MAX_DIFF, -m MAX_DIFF
                        Maximum sequence variations before querying sequencing accuracy
  -tmp PATH             Temporary folder location
  -output OUTFILE, -out OUTFILE, -o OUTFILE
                        Specify output file
  -threads THREADS, -t THREADS
                        Maximum threads
  -sort                 Run the "sortseq.py" script
  -all, -a              Test all sequences (VERY slow)
  -json                 Export results matrix as JSON file (Experimental)
```

#### Start your own analysis

With a `covid_genomes.fasta` as a multi-sequence fasta file (gzip compressed files are accepted) in the current directory.

```sh
docker run -i -t --rm -v "$(pwd)":"$(pwd)" -u $(id -u):$(id -g) covidmutants \
  -threads 15 \
  -forward TATgCCATTAgTgCAAAgAATAgAgCTCgCAC \
  -reverse GTAATTGGAAcAAGcAAATTcTATGGTGGTTG \
  -amplicon TATGCCATTAGTGCAAAGAATAGAGCTCGCACCGTAGCTGGTGTCTCTATCTGTAGTACTATGACCAATAGACAGTTTCATCAAAAATTATTGAAATCAATAGCCGCCACTAGAGGAGCTACTGTAGTAATTGGAACAAGCAAATTCTATGGTGGTTG \
  -o "$(pwd)"/output \
  -fasta "$(pwd)"/covid_genomes.fasta
```

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
...
```

#### \<output\>.fail.fa file

FAIL indicates that the primer sequences did not match the primers submitted. An extra file \<output\>.fail.fa is created with those sequences for manual validation.

#### \<output\>.issue.fa file

If the sequences (primer or amplicon) are found but have too many sequence variants (e.g. more than `-max`), the accuracy of their sequences is questioned.
The sequences are reported separately in the \<output\>.issue.fa file.

## Extra scripts

Three scripts are available.

* `typetempsort.py` extra script, parsing GISAID compressed exportation files into variant/month output FASTA files
* `oligomutk.py` the script run by the docker
* `sortseq.py` extra script, summarising `oligomutk.py` output.

#### `typetempsort.py` usage

```plaintext
usage: typetempsort.py [-h] -m METAPATH -f SEQUENCEPATH -s START_DATE -e END_DATE [-v VIRUS]

Parse GISAID bulk files into monthly parcels.

  -h, --help            show this help message and exit
  -i INPUTPATH, --input INPUTPATH
                        File path to the GISAID "Downloads > Genomic epidemilogy - nextregions" export file "hcov_xxxxxx_xxxx-xx-xx_xx-xx.tar.gz" file
  -s START_DATE, --start START_DATE
                        Start date (inclusive): eg, 2020-05 (for May 2020)
  -e END_DATE, --end END_DATE
                        End date (inclusive): eg, 2021-05 (for May 2021)
  -v VIRUS, --virus VIRUS
                        Virus subtype to be extracted: eg, Mu, Alpha, Beta, Gamma, Delta, GH/490R, Lambda, Omicron or ALL
```

###### Example

```sh
./typetempsort.py --input hcov_xxxxxx_xxxx-xx-xx_xx-xx.tar.gz \
  --start 2020-05 --end 2023-04 --virus Delta
```

#### `sortseq.py` usage

```plaintext
usage: sortseq.py [-h] -i INFILE [-o OUTFILE]

Summarise results.

  -h, --help            show this help message and exit
  -i INFILE, --in INFILE
                        input file (oligomutk.py output)
  -o OUTFILE, --out OUTFILE
                        result/summary file (TSV formatted)
```

###### Example

```sh
./sortseq.py -i docker_output.tsv -o output.summary.tsv
```


## Full pipeline (to reproduce the results in the publication)

#### Prepare the GIT repository and docker

```sh
git clone https://github.com/aicuramedical/covidmutants.git
cd covidmutants
docker build --rm=true -t covidmutants .
```

#### Collect the sequences

From GIDAID, manually download the metadata and sequences from the 2023-04-25 for oceania only: "Downloads" > "Genomic epidemilogy - nextregions" > "Oceania": `hcov_oceania_2023-04-25_08-29.tar.gz`. Then, split the file per month for the Omicron variant:

```sh
./typetempsort.py --input hcov_oceania_2023-04-25_08-29.tar.gz \
  --start 2021-12 --end 2023-04 --virus Omicron
```

Now you should have 17 files from `OMICRON_2021-12.fasta.gz` to `OMICRON_2023-04.fasta.gz`.

#### Run the full analysis

For each month/variant, run an analysis:

```sh
docker run -i -t --rm -v "$(pwd)":"$(pwd)" -u $(id -u):$(id -g) covidmutants \
  -threads 15 -trust \
  -forward TATgCCATTAgTgCAAAgAATAgAgCTCgCAC \
  -reverse GTAATTGGAAcAAGcAAATTcTATGGTGGTTG \
  -amplicon TATGCCATTAGTGCAAAGAATAGAGCTCGCACCGTAGCTGGTGTCTCTATCTGTAGTACTATGACCAATAGACAGTTTCATCAAAAATTATTGAAATCAATAGCCGCCACTAGAGGAGCTACTGTAGTAATTGGAACAAGCAAATTCTATGGTGGTTG \
  -o "$(pwd)"/output \
  -fasta "$(pwd)"/OMICRON_2022-11.fasta.gz
  -sort
```

The `output.tsv` file have the details results (produced by the `oligomutk.py` script), while the `output.summary.tsv` present the different mutations existing (produced by the `sortseq.py` script).

## Issues

If you have any problems with or questions about the scripts, please contact us through a [GitHub issue](https://github.com/aicuramedical/covidmutants/issues).
Any issue related to the scientific results themselves must be done directly with the authors.

## Contributing

You are invited to contribute new features, fixes, or updates, large or small; we are always thrilled to receive pull requests, and do our best to process them as fast as we can.

## License and distribution

The content of this project itself including the raw data and work are licensed under the [Creative Commons Attribution-ShareAlike 4.0 International License](http://creativecommons.org/licenses/by-sa/4.0/), and the source code presented is licensed under the [GPLv3 license](http://www.gnu.org/licenses/gpl-3.0.html).
