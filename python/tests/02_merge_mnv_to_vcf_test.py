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
Tests of the merge_mnv_to_vcf module
"""
import os

import pytest
from casmsmartphase import merge_mnv_to_vcf

TEST_INPUT = "test_data/test_input.vcf.gz"
EXP_OUTPUT_NORM = "test_data/test_exp_result.vcf"
TEST_OUTPUT = "test_data/test_output.vcf"
TEST_INPUT_FILTER = "test_data/test_input_filt_qual.vcf.gz"
EXP_OUTPUT_FILTER = "test_data/test_filt_qual_exp_result.vcf"
SMART_PHASE_OUTPUT = "test_data/sample.phased.output"
CUTOFF = 0.0
EXCLUDE = 2


def compare_files(file_a, file_b):
    """
    Utility method for comparison of files in txt format
    """
    if file_a == file_b:
        return True
    with open(file_a) as a:
        with open(file_b) as b:
            linea = a.readline()
            lineb = b.readline()
            assert linea == lineb
    return True


@pytest.mark.parametrize(
    "input,output,exp_output",
    [
        (TEST_INPUT, TEST_OUTPUT, EXP_OUTPUT_NORM),
        (TEST_INPUT_FILTER, TEST_OUTPUT, EXP_OUTPUT_FILTER),
    ],
)
def test_merge_mnv_to_vcf(input, output, exp_output):
    merge_mnv_to_vcf.run(input, output, SMART_PHASE_OUTPUT, CUTOFF, EXCLUDE)
    compare_files(output, exp_output)
    os.remove(output)
