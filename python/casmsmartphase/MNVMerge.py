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
#    1. The usage of a range of years within a copyright statement contained within
#    this distribution should be interpreted as being equivalent to a list of years
#    including the first and last year specified and all consecutive years between
#    them. For example, a copyright statement that reads ‘Copyright (c) 2005, 2007-
#    2009, 2011-2012’ should be interpreted as being identical to a statement that
#    reads ‘Copyright (c) 2005, 2007, 2008, 2009, 2011, 2012’ and a copyright
#    statement that reads ‘Copyright (c) 2005-2012’ should be interpreted as being
#    identical to a statement that reads ‘Copyright (c) 2005, 2006, 2007, 2008,
#    2009, 2010, 2011, 2012’."
#
#
"""
Python module for methods reading/parsing and merging MNVs in an
unflagged (not post processed) CaVEMan generated VCF
"""
import datetime
import os
import sys

import vcfpy

# Setup base variables for the VCF process line
base_vcf_process_key = "vcfProcessLog"
base_vcf_process_log = "<InputVCF=<{}>,InputVCFSource=<{}>,InputVCFParam=<{}>>"


def get_last_vcf_process_index(header_lines, key_prefix):
    """
    Find all header lines matching the prefix key_prefix.
    Where they do match, find the greatest index, so we can increment
    when generating a new header line for this process.
    Returns greatest index found in the header key.
    """
    matching_head_keys = set(
        line.key
        for line in (
            line
            for line in header_lines
            if isinstance(line, vcfpy.header.HeaderLine)
            and line.key.startswith(key_prefix)
        )
    )
    # Now find the maximum indices of the keys found
    indices = [
        int(key_part[1])
        for key_part in (key.rsplit(".", 1) for key in matching_head_keys)
        if len(key_part) == 2 and key_part[1].isdigit
    ]
    if indices:
        return max(indices)
    else:
        return None


class MNVMerge:
    """
    Class containing VCF parsing and MNV merging code
    """

    def __init__(self, vcfIn, vcfOut, run_script, arg_str):
        self.vcfinname = os.path.basename(vcfIn)
        self.vcfin = vcfpy.Reader.from_path(vcfIn)
        self.vcfout = vcfOut
        self.run_script = run_script
        self.arg_str = arg_str
        self.longest_MNV = 2

    def get_process_header_line(self, existing_head):
        """
        Generates a new vcfProvcess header line for this process.
        Uses the existing header to check whether we require an index
        and (UTC) date appended to the key
        """
        head_key = base_vcf_process_key
        index = get_last_vcf_process_index(existing_head.lines, head_key)
        # Test for key without date - gives all lines
        if index is not None and index > 0:
            now = datetime.datetime.utcnow()
            date_str = now.strftime("%Y%m%d")
            head_key += "_{}".format(date_str)
            # Test for key with current date to get precise index
            index = get_last_vcf_process_index(existing_head.lines, head_key)
            if index is None:
                index = 0
            head_key = head_key + "." + str(index + 1)
        new_process_line = vcfpy.HeaderLine(
            key=head_key,
            value=base_vcf_process_log.format(
                self.vcfinname, self.run_script, self.arg_str
            ),
        )
        return new_process_line

    def generate_new_increment_header(self, existing_line, n):
        """
        Taking a header line, and an int n, generates a copy of that header
        line with the key and description updates to include said
        incremental int
        """
        new_line = existing_line.copy()
        new_line.mapping["ID"] = new_line.mapping["ID"] + f"_{n}"
        new_line.mapping["Description"] = (
            new_line.mapping["Description"] + f" (MNV allele {n} in series)"
        )
        return new_line

    def parse_header_add_merge_and_process(self, writer_header):
        """
        Parse VCF header from input. Add a process line and MNV
        format lines
        """
        # Add a headerline to say this was refiltered with this tool
        process_head_line = self.get_process_header_line(writer_header)
        writer_header.add_line(process_head_line)
        info_lines = writer_header.get_lines("INFO")
        format_lines = writer_header.get_lines("FORMAT")
        # Add a X1..X2..Xn info and format header lines for
        lines_to_add = []
        for inf_line in info_lines:
            for i in range(1, self.longest_MNV + 1):
                lines_to_add.append(self.generate_new_increment_header(inf_line, i))
        for form_line in format_lines:
            for i in range(1, self.longest_MNV + 1):
                lines_to_add.append(self.generate_new_increment_header(form_line, i))

        for new_head_line in lines_to_add:
            writer_header.add_line(new_head_line)

        return writer_header

    def merge_snv_to_mnv(self, snv_list):
        """
        Merge snvs from list into a single variant and output to VCF
        """
        do_qual = 0
        do_filter = 0
        qual = 0
        filter = []
        if snv_list[0].QUAL:
            do_qual = 1
            print(
                "Found a QUAL value, output will be a mean of all QUAL.",
                file=sys.stderr,
                end="\n",
            )
        else:
            qual = snv_list[0].QUAL

        if snv_list[0].FILTER:
            do_filter = 1
            print(
                "Found a FILTER value, output will be all FILTERs encountered at all bases.",
                file=sys.stderr,
                end="\n",
            )
        else:
            filter = snv_list[0].FILTER

        chrom = snv_list[0].CHROM
        pos = snv_list[0].POS
        id = []  # list of SNV IDs?
        ref = ""
        alt_str = ""
        info = {}
        format = []
        # Setup new call objects
        calls_dict = {}

        for n, var in enumerate(snv_list, start=1):
            ref += str(var.REF)
            alt_str += str(var.ALT[0].value)

            # Add incremental format strings
            for f in var.FORMAT:
                format.append(f + f"_{n}")

            # Add incremented info
            for key, val in var.INFO.items():
                info[key + f"_{n}"] = val

            # Calls, should be a NORMAL and TUMOUR call and associated counts
            for call in var.calls:
                new_call = None

                # Check for existing new_call object
                data_to_add = dict()
                sample_nom = call.sample

                # Build new data from existing call data (add increment)
                for k, v in call.data.items():
                    k = k + f"_{n}"
                    data_to_add[k] = v
                if call.sample in calls_dict:
                    existing_call = calls_dict[call.sample]
                    exist_data = existing_call.data.copy()

                    # Append new data to previously processed
                    for k, v in data_to_add.items():
                        exist_data[k] = v
                    data_to_add = exist_data

                new_call = vcfpy.Call(sample_nom, data_to_add)
                calls_dict[call.sample] = new_call

            # If we want to append filters
            if do_filter:
                filter.append(var.FILTER)
            # If we want to append qualities to generate a mean

            if do_qual:
                qual += var.QUAL
            id.append(var.ID[0])
        alt = [vcfpy.Substitution(type_="MNV", value=alt_str)]

        # Make the quality a mean value for MNVs
        if do_qual:
            qual = qual / len(snv_list)

        # INFO    FORMAT  NORMAL  TUMOUR
        # DP=266;MP=1.0e+00;GP=4.6e-51;TG=GG/AGGGG;TP=1.0e+00;SG=GG/AAGGG;SP=9.4e-07      GT:FAZ:FCZ:FGZ:FTZ:RAZ:RCZ:RGZ:RTZ:PM   0|0:0:0:88:0:0:0:88:0:0.0e+00   0|1:6:0:40:0:4:0:40:0:1.1e-01
        # {'DP': 271, 'MP': 1.0, 'GP': 7.1e-49, 'TG': 'GG/AGGGG', 'TP': 1.0, 'SG': 'GG/AAGGG', 'SP': 1.6e-06}
        # {'DP': 266, 'MP': 1.0, 'GP': 4.6e-51, 'TG': 'GG/AGGGG', 'TP': 1.0, 'SG': 'GG/AAGGG', 'SP': 9.4e-07}
        # ['GT', 'FAZ', 'FCZ', 'FGZ', 'FTZ', 'RAZ', 'RCZ', 'RGZ', 'RTZ', 'PM']
        # ['GT', 'FAZ', 'FCZ', 'FGZ', 'FTZ', 'RAZ', 'RCZ', 'RGZ', 'RTZ', 'PM']
        # [Call('NORMAL', {'GT': '0|0', 'FAZ': 0, 'FCZ': 0, 'FGZ': 88, 'FTZ': 0, 'RAZ': 0, 'RCZ': 0, 'RGZ': 88, 'RTZ': 0, 'PM': 0.0}), Call('TUMOUR', {'GT': '0|1', 'FAZ': 6, 'FCZ': 0, 'FGZ': 40, 'FTZ': 0, 'RAZ': 4, 'RCZ': 0, 'RGZ': 40, 'RTZ': 0, 'PM': 0.11})]
        # [Call('NORMAL', {'GT': '0|0', 'FAZ': 0, 'FCZ': 1, 'FGZ': 88, 'FTZ': 0, 'RAZ': 1, 'RCZ': 0, 'RGZ': 88, 'RTZ': 0, 'PM': 0.0056}), Call('TUMOUR', {'GT': '0|1', 'FAZ': 6, 'FCZ': 0, 'FGZ': 41, 'FTZ': 0, 'RAZ': 5, 'RCZ': 1, 'RGZ': 40, 'RTZ': 0, 'PM': 0.12})]
        mnv = vcfpy.Record(
            chrom, pos, id, ref, alt, qual, filter, info, format, calls_dict.values()
        )

        if len(snv_list) > self.longest_MNV:
            self.longest_MNV = len(snv_list)
        return mnv

    def perform_mnv_merge(self):
        """
        Iterate through VCF records. Outputting a new VCF with
        requested filters removed.
        """
        reader = self.vcfin
        # Make a copy of the header
        writer_header = reader.header.copy()
        prev_snv = []
        prev_pos = 0
        prev_contig = ""
        vars_to_print = []
        for variant in reader:

            # If this variant is not adjacent to the previous
            if len(prev_snv) > 0 and (
                str(variant.CHROM) != prev_contig or int(variant.POS) > prev_pos + 1
            ):

                # Print any already adjacent SNVs as MNVs
                if len(prev_snv) > 1:
                    # MNVs require parsing and printing
                    mnv_rec = self.merge_snv_to_mnv(prev_snv)
                    vars_to_print.append(mnv_rec)
                    prev_snv.clear()

                else:
                    # length is only one so a single SNV to print back to VCF
                    vars_to_print.append(prev_snv.pop())

            prev_snv.append(variant)
            prev_contig = str(variant.CHROM)
            prev_pos = int(variant.POS)

        # Print any already adjacent SNVs as MNVs
        if len(prev_snv) > 1:
            # MNVs require parsing and printing
            mnv_rec = self.merge_snv_to_mnv(prev_snv)
            vars_to_print.append(mnv_rec)
            prev_snv.clear()

        else:
            # length is only one so a single SNV to print back to VCF
            vars_to_print.append(prev_snv.pop())

        writer_header = self.parse_header_add_merge_and_process(writer_header)

        writer = vcfpy.Writer.from_path(self.vcfout, writer_header)

        for vars in vars_to_print:
            writer.write_record(vars)

        writer.close()
