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
Merge adjacent SNVs in a CaVEMan called as MNVs by Smart-Phase into a
new VCF containing SNVs and MNVs
"""
import os

from casmsmartphase.MNVMerge import MNVMerge


def run(vcfin, output, smart_phased_output, cutoff, exclude, arg_str, bed=None):
    # Generate a merged VCF with possible MNVs
    # Open vcf reading module
    mnvmerge = MNVMerge(
        vcfin,
        output,
        smart_phased_output,
        cutoff,
        exclude,
        os.path.basename(__file__),
        arg_str,
        bed=None,
    )
    mnvmerge.perform_mnv_merge_to_vcf()
