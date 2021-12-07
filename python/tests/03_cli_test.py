# LICENSE
#
# Copyright (c) 2021
#
# Author: CASM/Cancer IT <cgphelp@sanger.ac.uk>
#
# This file is part of CASM-Smart-Phase.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
# 1. The usage of a range of years within a copyright statement contained within
# this distribution should be interpreted as being equivalent to a list of years
# including the first and last year specified and all consecutive years between
# them. For example, a copyright statement that reads ‘Copyright (c) 2005, 2007-
# 2009, 2011-2012’ should be interpreted as being identical to a statement that
# reads ‘Copyright (c) 2005, 2007, 2008, 2009, 2011, 2012’ and a copyright
# statement that reads ‘Copyright (c) 2005-2012’ should be interpreted as being
# identical to a statement that reads ‘Copyright (c) 2005, 2006, 2007, 2008,
# 2009, 2010, 2011, 2012’.
"""
Tests the cli module
"""
import sys

import pkg_resources  # part of setuptools
import pytest
from casmsmartphase.cli import cli
from click.testing import CliRunner

INPUT_VCF = "test_data/test_input.vcf.gz"
SPHASE_OUT = "test_data/sample.phased.output"
VERSION = pkg_resources.require("casmsmartphase")[0].version
EXP_BASE_HELP = """Usage: cli [OPTIONS] COMMAND [ARGS]...

Options:
  --version  Show the version and exit.
  --help     Show this message and exit.

Commands:
  generate-bed  Generate a bed file of adjacent SNVs in a VCF for...
  merge-mnvs    Merge MNVs parsed by smartphase into a CaVEMan SNV and MNV...
"""

EXP_GENERATE_BED_HELP = """Usage: cli generate-bed [OPTIONS]

  Generate a bed file of adjacent SNVs in a VCF for smartphase analysis

Options:
  --version                Show the version and exit.
  -f, --vcfin FILE         Path to input VCF file  [required]
  -o, --output output.bed  Path to write output bed file
  --markhz / --nomarkhz    Mark homozygous adjacent SNVs in the bed file output
                           (default - don't mark)
  --help                   Show this message and exit.
"""

EXP_MERGE_MNV_HELP = """Usage: cli merge-mnvs [OPTIONS]

  Merge MNVs parsed by smartphase into a CaVEMan SNV and MNV vcf file

Options:
  --version                       Show the version and exit.
  -f, --vcfin FILE                Path to input VCF file  [required]
  -o, --output output.vcf         Path to write output vcf file
  -p, --smart-phased-output sample.phased.output
                                  The phased output file from Smart-Phase
                                  [required]
  -c, --cutoff FLOAT              Exclude any MNVs with a phased score < cutoff
                                  [default: 0.0]
  -x, --exclude INTEGER           Exclude phased MNV if it matches any of the
                                  exclude flag bits
  -b, --bed FILE                  .bed file of regions used to run smartphase.
                                  If homozygous adjacent SNVs are marked in the
                                  file they will be output in the merged VCF as
                                  an MNV.
  --help                          Show this message and exit.
"""

runner = CliRunner()


def test_version():
    response = runner.invoke(cli, "--version")
    assert response.exit_code == 0
    assert response.output == f"cli, version {VERSION}\n"


def test_help():
    response = runner.invoke(cli, ["--help"])
    assert response.exit_code == 0
    assert response.output == EXP_BASE_HELP


def test_generate_bed():
    response = runner.invoke(cli, ["generate-bed", "--help"])
    assert response.output == EXP_GENERATE_BED_HELP
    assert response.exit_code == 0


def test_merge_mnvs():
    response = runner.invoke(cli, ["merge-mnvs", "--help"])
    assert response.output == EXP_MERGE_MNV_HELP
    assert response.exit_code == 0
