#!/usr/bin/env python
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
"""
merge_vcf_MNVs.py: Script to merge adjacent SNVs in a CaVEMan generated
VCF into a new VCF containing SNVs and merged MNVs in order to be
processed by Smart-phase
"""
# Core imports
import argparse
import os
import sys

import pkg_resources  # part of setuptools
from casmsmartphase.MNVMerge import MNVMerge

VERSION = pkg_resources.require("casmsmartphase")[0].version


# Define arguments and argument associated methods
def generate_arg_str(parsed_args):

    ag_str = ""
    idx = 0
    for arg in sorted(parsed_args.__dict__.keys()):
        val = parsed_args.__dict__[arg]
        if not val:
            continue
        if arg == "vcfin" or arg == "vcfout" or arg == "bedfile":
            val = os.path.basename(val)

        if arg == "flagremove":
            inner_idx = idx
            for v in val:
                if inner_idx > 0:
                    ag_str += ","
                ag_str += "{}={}".format(arg, v)
                inner_idx += 1
        else:
            if idx > 0:
                ag_str += ","
            ag_str += "{}={}".format(arg, val)
        idx += 1
    return ag_str


parser = argparse.ArgumentParser(
    description="""
                  Merge adjacent SNVs in a VCF into
                  a new VCF containing SNVs and merged MNVs
                  """,
    prog="merge_vcf_MNVs.py",
)

parser.add_argument(
    "-f",
    "--vcfin",
    dest="vcfin",
    metavar="input.vcf[.gz]",
    help="Path to input VCF file",
    required=True,
)

parser.add_argument(
    "-p",
    "--smart-phased-output",
    dest="spout",
    help="The phased output file from Smart-Phase",
    required=True,
    metavar="sample.phased.output",
)

parser.add_argument(
    "-o",
    "--output",
    metavar="output.vcf",
    help="path to write output VCF file",
    required=False,
    default="<vcfin>.MNV.vcf",
)

parser.add_argument(
    "-c",
    "--cutoff-score",
    default=0,
    type=float,
    dest="cutoff",
    help="Exclude any MNVs with a phased score < cutoff",
    required=False,
)

parser.add_argument(
    "-x",
    "--exclude-flag",
    default=2,  # Default to exclude 'trans' phased MNVs
    type=int,
    dest="exclude",
    help="Exclude phased MNV if it matches any of the exclude flag bits",
    required=False,
)

parser.add_argument(
    "-v",
    "--version",
    dest="version",
    action="version",
    help="Print version information",
    version="%(prog)s " + VERSION,
)

args = parser.parse_args()
arg_str = generate_arg_str(args)

# Generate a merged VCF with possible MNVs
# Open vcf reading module
mnvmerge = MNVMerge(
    args.vcfin,
    args.output,
    args.spout,
    args.cutoff,
    args.exclude,
    os.path.basename(__file__),
    arg_str,
)
mnvmerge.perform_mnv_merge_to_vcf()
