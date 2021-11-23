# CASM-Smart-Phase/python

- [Installation](#installation)
  - [Without Docker image](#without-docker-image)
- [Python Utility Commands](#python-utility-commands)
  - [generate-bed](#generate-bed)
  - [merge-mnvs](#merge-mnvs)

## Installation

### Without Docker image

- Requires >= python 3.7

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

<!-- markdownlint-disable MD036 -->

**Alternatively use the CASM-Smart-Phase Docker container**

## Python Utility Commands

### generate-bed

Generate a bed file of adjacent SNVs in a VCF file. Allows for targeted analysis in smart-phase using the `-g` option

```bash
$ casmsmartphase generate-bed --help
Usage: casmsmartphase generate-bed [OPTIONS]

  Generate a bed file of adjacent SNVs in a VCF for smartphase analysis

Options:
  --version                Show the version and exit.
  -f, --vcfin FILE         Path to input VCF file  [required]
  -o, --output output.bed  Path to write output bed file
  --markhz / --nomarkhz    Mark homozygous adjacent SNVs in the bed file
                           output (default - don't mark)
  --help                   Show this message and exit.
```

### merge-mnvs

Merge MNVs from original VCF using the output from smart-phase

```bash
$ casmsmartphase merge-mnvs --help
Usage: casmsmartphase merge-mnvs [OPTIONS]

  Merge MNVs parsed by smartphase into a CaVEMan SNV and MNV vcf file

Options:
  -f, --vcfin FILE                Path to input VCF file  [required]
  -o, --output output.vcf         Path to write output vcf file
  -p, --smart-phased-output sample.phased.output
                                  The phased output file from Smart-Phase
                                  [required]
  -c, --cutoff FLOAT              Exclude any MNVs with a phased score <
                                  cutoff
  -x, --exclude INTEGER           Exclude phased MNV if it matches any of the
                                  exclude flag bits
  --help                          Show this message and exit.
```
