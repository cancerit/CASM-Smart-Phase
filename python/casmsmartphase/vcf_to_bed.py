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
Merge adjacent SNVs in a CaVEMan generated
VCF into a new VCF containing SNVs and merged MNVs in order to be
processed by Smart-phase
"""
import vcfpy


def parse_vcf(reader, outfile, markhz=False):
    print(f"markhz {markhz}")
    prev_contig = ""
    prev_snv = []
    all_is_het = -1  # -1 unset, 0 hom, 1 het
    prev_pos = 0
    for variant in reader:

        # If this variant is not adjacent to the previous
        if len(prev_snv) > 0 and (
            variant.CHROM != prev_contig or int(variant.POS) > prev_pos + 1
        ):
            # Print any already adjacent SNVs as MNVs
            if len(prev_snv) > 1:
                # MNVs print possible MNV location to bed file
                bed_str = (
                    f"{prev_snv[0].CHROM}\t{prev_snv[0].POS-1}\t{prev_snv[-1].POS}"
                )
                if all_is_het == 0 and markhz:
                    bed_str = bed_str + f"\thom"
                print(
                    bed_str,
                    file=outfile,
                )
            prev_snv.clear()
            all_is_het = -1

        prev_snv.append(variant)
        calculated_prev_het = 0
        if (variant.call_for_sample["TUMOUR"]).is_het or all_is_het == 1:
            calculated_prev_het = 1
        all_is_het = calculated_prev_het
        prev_contig = str(variant.CHROM)
        prev_pos = int(variant.POS)

    # Print any already adjacent SNVs as MNVs
    if len(prev_snv) > 1:
        # MNVs pront possible MNV location to bed file
        bed_str = f"{prev_snv[0].CHROM}\t{prev_snv[0].POS-1}\t{prev_snv[-1].POS}"
        if all_is_het == 0 and markhz:
            bed_str = bed_str + f"\thom"
        print(
            bed_str,
            file=outfile,
        )
    prev_snv.clear()


def run(vcfin, output, markhz):
    # Run through input VCF file and output any bed locations
    """
    Iterate through VCF records. Outputting a new VCF with
    requested filters removed.
    """
    reader = vcfpy.Reader.from_path(vcfin)
    with open(output, "w") as outfile:
        parse_vcf(reader, outfile, markhz)
