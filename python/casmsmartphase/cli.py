##  LICENSE
# Copyright (c) 2021
# Author: CASM/Cancer IT <cgphelp@sanger.ac.uk>
# This file is part of CASM-Smart-Phase.
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
# 1. The usage of a range of years within a copyright statement contained within
# this distribution should be interpreted as being equivalent to a list of years
# including the first and last year specified and all consecutive years between
# them. For example, a copyright statement that reads ‘Copyright (c) 2005, 2007-
# 2009, 2011-2012’ should be interpreted as being identical to a statement that
# reads ‘Copyright (c) 2005, 2007, 2008, 2009, 2011, 2012’ and a copyright
# statement that reads ‘Copyright (c) 2005-2012’ should be interpreted as being
# identical to a statement that reads ‘Copyright (c) 2005, 2006, 2007, 2008,
# 2009, 2010, 2011, 2012’.
from functools import wraps

import click
import pkg_resources  # part of setuptools
from casmsmartphase import merge_mnv_to_vcf
from casmsmartphase import vcf_to_bed

HELP_VCF_IN = "Path to input VCF file"
HELP_EXCLUDE = "Exclude phased MNV if it matches any of the exclude flag bits"
HELP_CUTOFF = "Exclude any MNVs with a phased score < cutoff"
HELP_OUTPUT_BED = "Path to write output bed file"
HELP_OUTPUT_VCF = "Path to write output vcf file"
HELP_SPHASE_OUT = "The phased output file from Smart-Phase"


def _file_exists():
    return click.Path(
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
        resolve_path=True,
    )


def common_params(f):
    @click.version_option(pkg_resources.require(__name__.split(".")[0])[0].version)
    @click.option(
        "-f",
        "--vcfin",
        required=True,
        default=None,
        type=_file_exists(),
        help=HELP_VCF_IN,
    )
    @wraps(f)
    def wrapper(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper


@click.group()
@click.version_option(pkg_resources.require(__name__.split(".")[0])[0].version)
def cli():
    pass


@cli.command()
@common_params
@click.option(
    "-o",
    "--output",
    metavar="output.bed",
    help=HELP_OUTPUT_BED,
    required=False,
    default="mnv.bed",
)
def generate_bed(*args, **kwargs):
    """
    Generate a bed file of adjacent SNVs in a VCF for smartphase analysis
    """
    vcf_to_bed.run(*args, **kwargs)


@cli.command()
@common_params
@click.option(
    "-o",
    "--output",
    metavar="output.vcf",
    help=HELP_OUTPUT_VCF,
    required=False,
    default="<vcfin>.MNV.vcf",
)
@click.option(
    "-p",
    "--smart-phased-output",
    help=HELP_SPHASE_OUT,
    required=True,
    metavar="sample.phased.output",
)
@click.option(
    "-c",
    "--cutoff",
    default=0.0,
    type=float,
    help=HELP_CUTOFF,
    required=False,
)
@click.option(
    "-x",
    "--exclude",
    default=2,  # Default to exclude 'trans' phased MNVs
    type=int,
    help=HELP_EXCLUDE,
    required=False,
)
def merge_mnvs(*args, **kwargs):
    """
    Merge MNVs parsed by smartphase into a CaVEMan SNV and MNV vcf file
    """
    merge_mnv_to_vcf.run(*args, **kwargs)
