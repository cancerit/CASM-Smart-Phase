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
Tests of the MNVMerge module
"""
import os

import pytest
import vcfpy
from casmsmartphase.MNVMerge import get_last_vcf_process_index
from casmsmartphase.MNVMerge import MNVMerge
from casmsmartphase.MNVMerge import parse_homs_bed_to_dict
from casmsmartphase.MNVMerge import parse_sphase_output

INPUT_VCF = "test_data/test_input.vcf.gz"
TRINUC_INPUT_VCF = "test_data/test_input_trinuc.vcf.gz"
FILT_QUAL_INPUT_VCF = "test_data/test_input_filt_qual.vcf.gz"
OUTPUT_VCF = "test_data/test_output.vcf"
FILT_QUAL_EXP_RES_VCF = "test_data/test_filt_qual_exp_result.vcf"
EXP_RES_VCF = "test_data/test_exp_result.vcf"
TRINUC_EXP_RES_VCF = "test_data/test_exp_result_trinuc.vcf"
RUN_SCRIPT = "pytest_MNVMerge"
ARG_STR = "x=test_Arg_str"
SPOUT = "test_data/sample.phased.output"
SPOUT_TRINUC = "test_data/sample.phased.trinuc.output"
BAD_SPOUT = "test_data/bad_sample.phased.output"
SPOUT_EXCEPT = "test_data/sample.phased.except.output"
BED_INPUT_HOM = "test_data/expected_output_hethom.bed"
BED_INPUT_NOHOM = "test_data/expected_output.bed"
SPOUT_TRINUC_2 = "test_data/test_phase_triplet.out"
CUTOFF = 0.0
EXCLUDE = 2


def compare_variants(vcf1, vcf2):
    """
    Utility method for comparison of VCF variants in vcfpy format
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
    Utility method for comparison of VCF headers in the vcfpy format
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
    "bed_file, exp_result",
    [
        (BED_INPUT_NOHOM, {}),
        (
            BED_INPUT_HOM,
            {
                "chr1": [
                    (1866692, 1866693, True),
                ],
                "chr3": [
                    (45636146, 45636147, True),
                ],
            },
        ),
    ],
)
def test_parse_homs_bed_to_dict(bed_file, exp_result):
    assert parse_homs_bed_to_dict(bed_file) == exp_result


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
    "sphaseout,cutoff,exclude_flags,hom_dict,exp_result",
    [
        (SPOUT, 0.0, 2, {}, ({"chr1": {1627262: 1627263}}, 2)),
        (SPOUT, 0.1, 1, {}, ({}, 1)),
        (SPOUT_EXCEPT, 0.0, 2, {}, ({"chr1": {1627262: 1627263}}, 2)),
        (SPOUT_TRINUC, 0.0, 2, {}, ({"chr12": {9420710: 9420713}}, 4)),
        (
            SPOUT_TRINUC,
            0.0,
            2,
            {"chr1": [(1627262, 1627263, "hom")]},
            ({"chr1": {1627262: 1627263}, "chr12": {9420710: 9420713}}, 4),
        ),
        (
            SPOUT_TRINUC,
            0.0,
            2,
            {"chr1": [(1627262, 1627269, "hom")]},
            ({"chr1": {1627262: 1627269}, "chr12": {9420710: 9420713}}, 8),
        ),
        (
            SPOUT_TRINUC_2,
            0.0,
            2,
            {},
            (
                {
                    "chr17": {42760364: 42760366},
                },
                3,
            ),
        ),
    ],
)
def test_parse_sphase_output(sphaseout, cutoff, exclude_flags, hom_dict, exp_result):
    assert parse_sphase_output(sphaseout, cutoff, exclude_flags, hom_dict) == exp_result


def test_parse_sphase_output_err():
    parse_sphase_output(BAD_SPOUT, CUTOFF, EXCLUDE, {})


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
    merge_obj = MNVMerge(
        INPUT_VCF, OUTPUT_VCF, SPOUT, CUTOFF, EXCLUDE, RUN_SCRIPT, ARG_STR
    )
    new_header_line = merge_obj.generate_new_increment_header(existing_line, n)
    assert new_header_line == exp_line


@pytest.mark.parametrize(
    "invcf,exp_res,spout",
    [
        (INPUT_VCF, EXP_RES_VCF, SPOUT),
        (FILT_QUAL_INPUT_VCF, FILT_QUAL_EXP_RES_VCF, SPOUT),
        (TRINUC_INPUT_VCF, TRINUC_EXP_RES_VCF, SPOUT_TRINUC),
    ],
)
def test_perform_mnv_merge(invcf, exp_res, spout):
    merge_obj = MNVMerge(invcf, OUTPUT_VCF, spout, CUTOFF, EXCLUDE, RUN_SCRIPT, ARG_STR)
    merge_obj.perform_mnv_merge_to_vcf()
    assert compare_vcf_files(OUTPUT_VCF, exp_res)
    os.remove(OUTPUT_VCF)
