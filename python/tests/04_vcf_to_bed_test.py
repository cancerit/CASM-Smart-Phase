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
Tests of the vcf_to_bed module
"""
import os
import sys

import pytest
from casmsmartphase import vcf_to_bed

TEST_INPUT = "test_data/test_input.vcf.gz"
EXP_OUTPUT = "test_data/expected_output.bed"
TEST_OUTPUT = "test_data/test_output.bed"
TEST_INPUT_HOM = "test_data/test_input_hethom.vcf.gz"
EXP_OUTPUT_HOM = "test_data/expected_output_hethom.bed"


def compare_files(file_a, file_b):
    """
    Utility method for comparison of files in txt format
    """
    if file_a == file_b:
        return True
    print(file_a, file_b, file=sys.stderr)
    with open(file_a) as a:
        with open(file_b) as b:
            linesa = a.readlines()
            linesb = b.readlines()
            assert linesa == linesb
    return True


@pytest.mark.parametrize(
    "input,output,markhom,exp_out",
    [
        (
            TEST_INPUT,
            TEST_OUTPUT,
            False,
            EXP_OUTPUT,
        ),
        (
            TEST_INPUT_HOM,
            TEST_OUTPUT,
            True,
            EXP_OUTPUT_HOM,
        ),
    ],
)
def test_vcf_to_bed_run(input, output, markhom, exp_out):
    vcf_to_bed.run_parse(input, output, markhom)
    compare_files(exp_out, output)
    os.remove(output)
