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
import os
from functools import wraps

import click
import pkg_resources  # part of setuptools
from casmsmartphase import merge_mnv_to_vcf
from casmsmartphase import vcf_to_bed

CUTOFF_DEFAULT = 0.0
HELP_VCF_IN = "Path to input VCF file"
HELP_EXCLUDE = "Exclude phased MNV if it matches any of the exclude flag bits"
HELP_CUTOFF = (
    f"Exclude any MNVs with a phased score < cutoff [default: {CUTOFF_DEFAULT}]"
)
HELP_OUTPUT_BED = "Path to write output bed file"
HELP_OUTPUT_HZ_BED = (
    "Mark homozygous adjacent SNVs in the bed file output (default - don't mark)"
)
HELP_OUTPUT_VCF = "Path to write output vcf file"
HELP_SPHASE_OUT = "The phased output file from Smart-Phase"
HELP_BED_REGIONS = """.bed file of regions used to run smartphase.
                    If homozygous adjacent SNVs are marked in the file they will be output in the merged VCF as an MNV."""
FILEPATH_INPUTS = ["vcfin", "output", "smart_phased_output"]


def _file_exists():
    return click.Path(
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
        resolve_path=True,
    )


def generate_arg_string(*args, **kwargs):
    ag_str = ""
    idx = 0
    for key, item in kwargs.items():
        if key in FILEPATH_INPUTS:
            item = os.path.basename(item)
        if idx > 0:
            ag_str += ","
        ag_str += f"{key}={item}"
        idx += 1
    return ag_str


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
)
@click.option("--markhz/--nomarkhz", help=HELP_OUTPUT_HZ_BED, default=False)
def generate_bed(*args, **kwargs):
    """
    Generate a bed file of adjacent SNVs in a VCF for smartphase analysis
    """
    vcf_to_bed.run_parse(*args, **kwargs)


@cli.command()
@common_params
@click.option(
    "-o",
    "--output",
    metavar="output.vcf",
    help=HELP_OUTPUT_VCF,
    required=False,
    default="output.MNV.vcf",
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
    default=CUTOFF_DEFAULT,
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
@click.option(
    "-b",
    "--bed",
    required=False,
    default=None,
    type=_file_exists(),
    help=HELP_BED_REGIONS,
)
def merge_mnvs(*args, **kwargs):
    """
    Merge MNVs parsed by smartphase into a CaVEMan SNV and MNV vcf file
    """
    arg_str = generate_arg_string(*args, **kwargs)
    kwargs["arg_str"] = arg_str
    merge_mnv_to_vcf.run(*args, **kwargs)
