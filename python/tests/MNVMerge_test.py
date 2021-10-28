#   LICENSE
# Copyright (c) 2018-2019 Genome Research Ltd.
# Author: Cancer Genome Project cgphelp@sanger.ac.uk
#
#
# This file is part of vcf_flag_modifier.
#
# vcf_flag_modifier is free software: you can redistribute it and/or modify it under
# the terms of the GNU Affero General Public License as published by the Free
# Software Foundation; either version 3 of the License, or (at your option) any
# later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more
# details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#
#  1. The usage of a range of years within a copyright statement contained within
#  this distribution should be interpreted as being equivalent to a list of years
#  including the first and last year specified and all consecutive years between
#  them. For example, a copyright statement that reads ‘Copyright (c) 2005, 2007-
#  2009, 2011-2012’ should be interpreted as being identical to a statement that
#  reads ‘Copyright (c) 2005, 2007, 2008, 2009, 2011, 2012’ and a copyright
#  statement that reads ‘Copyright (c) 2005-2012’ should be interpreted as being
#  identical to a statement that reads ‘Copyright (c) 2005, 2006, 2007, 2008,
#  2009, 2010, 2011, 2012’."
#
#
"""
Tests of the MNVMerge module
"""
import os

import pytest
import vcfpy
from casmsmartphase.MNVMerge import get_last_vcf_process_index
from casmsmartphase.MNVMerge import MNVMerge

input_vcf = "test_data/test_input.vcf.gz"
output_vcf = "test_data/test_output.vcf"
exp_res_vcf = "test_data/test_exp_result.vcf"
run_script = "pytest_MNVMerge"
arg_str = "x=test_Arg_str"


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


def generate_new_increment_header(existing_line, n):
    assert 1 == 0

    # def generate_new_increment_header(self, existing_line, n):
    #     """
    #     Taking a header line, and an int n, generates a copy of that header
    #     line with the key and description updates to include said
    #     incremental int
    #     """
    #     new_line = existing_line.copy()
    #     new_line.mapping["ID"] = new_line.mapping["ID"] + f"_{n}"
    #     new_line.mapping["Description"] = (
    #         new_line.mapping["Description"] + f" (MNV allele {n} in series)"
    #     )
    #     return new_line


def test_perform_mnv_merge():
    merge_obj = MNVMerge(input_vcf, output_vcf, run_script, arg_str)
    merge_obj.perform_mnv_merge()
    assert compare_vcf_files(output_vcf, exp_res_vcf)
    os.remove(output_vcf)
