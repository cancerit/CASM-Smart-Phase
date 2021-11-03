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
Tests of the MNVMerge module
"""
import os

import pytest
import vcfpy
from casmsmartphase.MNVMerge import get_last_vcf_process_index
from casmsmartphase.MNVMerge import MNVMerge

INPUT_VCF = "test_data/test_input.vcf.gz"
OUTPUT_VCF = "test_data/test_output.vcf"
EXP_RES_VCF = "test_data/test_exp_result.vcf"
RUN_SCRIPT = "pytest_MNVMerge"
ARG_STR = "x=test_Arg_str"


def compare_variants(vcf1, vcf2):
    """
    Utility method for comparison of VCF variants in the cyvcf2 format
    """
    for var1, var2 in zip(vcf1, vcf2):
        assert var1.CHROM == var2.CHROM
        assert var1.POS == var2.POS
        assert var1.REF == var2.REF
        assert var1.ALT == var2.ALT
        assert var1.ID == var2.ID
        assert var1.QUAL == var2.QUAL
        assert var1.FILTER == var2.FILTER
        assert var1.FORMAT == var2.FORMAT
        assert dict(var1.INFO) == dict(var2.INFO)


def compare_vcf_header(vcf1, vcf2):
    """
    Utility method for comparison of VCF headers in the cyvcf2 format
    """
    # assert vcf1.header.lines == vcf2.header.lines

    assert list(vcf1.header.get_lines("contig")) == list(
        vcf2.header.get_lines("contig")
    )
    assert list(vcf1.header.get_lines("INFO")) == list(vcf2.header.get_lines("INFO"))
    assert list(vcf1.header.get_lines("FORMAT")) == list(
        vcf2.header.get_lines("FORMAT")
    )
    assert list(vcf1.header.get_lines("FILTER")) == list(
        vcf2.header.get_lines("FILTER")
    )


def compare_vcf_files(vcf_a, vcf_b):
    """
    Utility method for comparison of VCF files in the vcfpy format
    """
    if vcf_a == vcf_b:
        return True
    vcf1 = vcfpy.Reader.from_path(vcf_a)
    vcf2 = vcfpy.Reader.from_path(vcf_b)
    compare_vcf_header(vcf1, vcf2)
    compare_variants(vcf1, vcf2)
    return True


@pytest.mark.parametrize(
    "in_head,key_prefix,exp_idx",
    [
        ([vcfpy.HeaderLine("vcfProcessLog", "BLAH")], "vcfProcessLog", None),
        (
            [
                vcfpy.HeaderLine("vcfProcessLog_20160212.1", "BLAH1"),
                vcfpy.HeaderLine("vcfProcessLog_20160212.2", "BLAH2"),
            ],
            "vcfProcessLog",
            2,
        ),
        (
            [
                vcfpy.HeaderLine("vcfProcessLog_20151218.1", "BLAH1"),
                vcfpy.HeaderLine("vcfProcessLog_20160212.1", "BLAH2"),
            ],
            "vcfProcessLog",
            1,
        ),
        ([], "vcfProcessLog", None),
    ],
)
def test_get_last_vcf_process_index(in_head, key_prefix, exp_idx):
    assert get_last_vcf_process_index(in_head, key_prefix) == exp_idx


@pytest.mark.parametrize(
    "existing_line, n, exp_line",
    [
        (
            vcfpy.FormatHeaderLine.from_mapping(
                {
                    "ID": "TEST",
                    "Number": 1,
                    "Type": "Integer",
                    "Description": "Test_Description",
                }
            ),
            1,
            vcfpy.FormatHeaderLine.from_mapping(
                {
                    "ID": "TEST_1",
                    "Number": 1,
                    "Type": "Integer",
                    "Description": "Test_Description (MNV allele 1 in series)",
                }
            ),
        ),
        (
            vcfpy.InfoHeaderLine.from_mapping(
                {
                    "ID": "TEST",
                    "Number": 1,
                    "Type": "String",
                    "Description": "Test_Description",
                }
            ),
            2,
            vcfpy.InfoHeaderLine.from_mapping(
                {
                    "ID": "TEST_2",
                    "Number": 1,
                    "Type": "String",
                    "Description": "Test_Description (MNV allele 2 in series)",
                }
            ),
        ),
    ],
)
def test_generate_new_increment_header(existing_line, n, exp_line):
    merge_obj = MNVMerge(INPUT_VCF, OUTPUT_VCF, RUN_SCRIPT, ARG_STR)
    new_header_line = merge_obj.generate_new_increment_header(existing_line, n)
    assert new_header_line == exp_line


def test_perform_mnv_merge():
    merge_obj = MNVMerge(INPUT_VCF, OUTPUT_VCF, RUN_SCRIPT, ARG_STR)
    merge_obj.perform_mnv_merge()
    assert compare_vcf_files(OUTPUT_VCF, EXP_RES_VCF)
    os.remove(OUTPUT_VCF)
