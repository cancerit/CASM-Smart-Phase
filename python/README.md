# CASM-Smart-Phase/python

- [Installation](#installation)
  - [Without Docker image](#without-docker-image)
- [Python Utility Scripts](#python-utility-scripts)
  - [scripts/generate_MNV.bed.py](#scriptsgenerate_mnvbedpy)
  - [scripts/merge_vcf_MNVs.py](#scriptsmerge_vcf_mnvspy)

## Installation

### Without Docker image

- Requires >= python 3.6

```bash
# Setup a python virtual environment.
python -m venv /path/to/venv
# Activate python environment
source /path/to/venv/bin/activate
# CD into the python directory
cd CASM-Smart-Phase/python
# Install CASM-Smart-Phase python
pip install ./
```

**Alternatively use the Docker container**

## Python Utility Scripts

### [scripts/generate_MNV.bed.py](scripts/generate_MNV_bed.py)

Generate a bed file of adjacent SNVs in a VCF file. Allows for targeted analysis in smart-phase idsing the `-g` option

```bash
$ generate_MNV_bed.py -h
usage: generate_MNV_bed.py [-h] [-v] -f input.vcf[.gz] [-o output.bed]

Generate a bed file of adjacent SNVs in a VCF for smartphase analysis

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         Print version information
  -f input.vcf[.gz], --vcfin input.vcf[.gz]
                        Path to input VCF file
  -o output.bed, --output output.bed
```

### [scripts/merge_vcf_MNVs.py](scripts/merge_vcf_MNVs.py)

Merge MNVs from original VCF using the output from smart-phase

```bash
$ merge_vcf_MNVs.py -h
usage: merge_vcf_MNVs.py [-h] -f input.vcf[.gz] -p sample.phased.output
                         [-o output.vcf] [-c CUTOFF] [-x EXCLUDE] [-v]

Merge adjacent SNVs in a VCF into a new VCF containing SNVs and merged MNVs

optional arguments:
  -h, --help            show this help message and exit
  -f input.vcf[.gz], --vcfin input.vcf[.gz]
                        Path to input VCF file
  -p sample.phased.output, --smart-phased-output sample.phased.output
                        The phased output file from Smart-Phase
  -o output.vcf, --output output.vcf
                        path to write output VCF file
  -c CUTOFF, --cutoff-score CUTOFF
                        Exclude any MNVs with a phased score < cutoff
  -x EXCLUDE, --exclude-flag EXCLUDE
                        Exclude phased MNV if it matches any of the exclude
                        flag bits
  -v, --version         Print version information
```
