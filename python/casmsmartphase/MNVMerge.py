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
Python module for methods reading/parsing and merging MNVs in an
unflagged (not post processed) CaVEMan generated VCF
"""
import datetime
import logging
import os
import re
from itertools import groupby
from typing import Dict
from typing import List
from typing import Optional
from typing import Tuple

LOGGER = logging.getLogger(__name__)

import vcfpy

# Setup base variables for the VCF process line
BASE_VCF_PROCESS_KEY = "vcfProcessLog"
BASE_VCF_PROCESS_LOG = "<InputVCF=<{}>,InputVCFSource=<{}>,InputVCFParam=<{}>>"


def parse_homs_bed_to_dict(bed_file: str) -> Dict:
    bed_entries_by_contig = dict()
    with open(bed_file, "r") as read_bed:
        line = read_bed.readline()
        while line:
            line = line.rstrip()
            split_line = line.split("\t")
            hom = False if len(split_line) == 3 else True
            if hom:
                bed_entry = (int(split_line[1]) + 1, int(split_line[2]), hom)
                if split_line[0] not in bed_entries_by_contig:
                    bed_entries_by_contig[split_line[0]] = []
                bed_entries_by_contig[split_line[0]].append(bed_entry)
            line = read_bed.readline()
    return bed_entries_by_contig


def get_last_vcf_process_index(
    header_lines: List[vcfpy.header.HeaderLine], key_prefix: str
) -> Optional[int]:
    """
    Find all header lines matching the prefix key_prefix.
    Where they do match, find the greatest index, so we can increment
    when generating a new header line for this process.
    Returns greatest index found in the header key.
    """
    matching_head_keys = set(
        [
            line.key
            for line in header_lines
            if isinstance(line, vcfpy.header.HeaderLine)
            and line.key.startswith(key_prefix)
        ]
    )
    # Now find the maximum indices of the keys found
    max_index = None
    for key_part in (key.rsplit(".", 1) for key in matching_head_keys):
        if (
            len(key_part) == 2
            and key_part[1].isdigit
            and (max_index is None or int(key_part[1]) > max_index)
        ):
            max_index = int(key_part[1])
    return max_index


def all_equal(iterable):
    g = groupby(iterable)
    return next(g, True) and not next(g, False)


def parse_sphase_output(
    sphaseout: str, cutoff: float, exclude_flags: int, hom_bed_parsed: Dict
) -> Tuple[Dict, int]:
    mnvs = {}
    max_len = 1
    with open(sphaseout, "r") as readspout:
        while True:
            line = readspout.readline()
            if not line or line.startswith("Denovo count"):
                break
            line = line.rstrip()
            try:
                (mnv_id, pair1, pair2, flag, confidence) = re.split(r"\s+", line, 5)
                # (_id_contig, id_start_region, _id_stop) = mnv_id.split("-")
                # print(_id_contig, id_start_region, _id_stop)
                if float(confidence) < cutoff or int(flag) & exclude_flags:
                    continue
                (contig, startpos, _tmp) = pair1.split("-", maxsplit=2)
                (contig, endpos, _tmp) = pair2.split("-", maxsplit=2)
                startpos = int(startpos)
                endpos = int(endpos)
                if startpos + 1 != endpos:
                    # Skip as non-adjacent pair test
                    continue
                if not contig in mnvs:
                    mnvs[contig] = {}
                # Check for adjacent MNV
                if startpos in mnvs[contig].values():
                    # Check if startpos in mnv_id is already stored and see if these are adjacent to an already recorded MNV
                    # Find key for mnv that adjoins this one
                    key = None
                    for k, (end_pos, conf_score_list) in mnvs[contig].items():
                        if end_pos == startpos:
                            key = k
                            break
                    # Check if current end_pos is adjacent, if so, extend this MNV
                    del (mnvs[contig])[key]
                    mnvs[contig][key] = (endpos, conf_score_list.append(confidence))
                    print(mnvs)
                    mnv_len = (endpos - key) + 1
                    if max_len < mnv_len:
                        max_len = mnv_len
                else:
                    # Otherwise this is a new MNV
                    mnvs[contig][startpos] = (endpos, [confidence])
                    mnv_len = (endpos - startpos) + 1
                    if max_len < mnv_len:
                        max_len = mnv_len

            # Possibly a non phased entry, check for length 2 when split before erroring
            except ValueError as err:  # Possibly a non phased entry, check for length 2 when split before erroring
                if len(re.split(r"\s+", line)) != 2:
                    raise ValueError(
                        f"Error encountered parsing smart-phase output at line {line}.\nOriginal error {err}"
                    )
                LOGGER.info(
                    f"Skipping line of only 2 items, not a phased variant {line}"
                )
                continue
    # Add the mnv's that are hom to the MNV list
    if hom_bed_parsed:
        for contig in hom_bed_parsed.keys():
            if contig not in mnvs:
                mnvs[contig] = {}
            for store in hom_bed_parsed[contig]:
                # (start, stop, hom)
                (start, stop, _hom) = store
                mnv_len = (stop - start) + 1
                mnvs[contig][start] = (stop, ["1.0"] * mnv_len)
                if max_len < mnv_len:
                    max_len = mnv_len
            mnvs[contig] = {
                k: mnvs[contig][k] for k in sorted(mnvs[contig])
            }  # Sort by position
        mnvs = {k: mnvs[k] for k in sorted(mnvs)}  # Sort by contig
    return mnvs, max_len


class MNVMerge:
    """
    Class containing VCF parsing and MNV merging code
    """

    def __init__(
        self,
        vcfIn: str,
        vcfOut: str,
        spout: str,
        cutoff: float,
        exclude: int,
        run_script: str,
        arg_str: str,
        bed=None,
    ):
        self.vcfinname = os.path.basename(vcfIn)
        self.vcfin = vcfpy.Reader.from_path(vcfIn)
        self.vcfout = vcfOut
        self.spout = spout
        self.cutoff = cutoff
        self.run_script = run_script
        self.exclude_flags = exclude
        self.arg_str = arg_str
        self.longest_MNV = 2
        self.bed = bed

    def get_process_header_line(
        self, existing_head: vcfpy.Header
    ) -> vcfpy.header.HeaderLine:
        """
        Generates a new vcfProvcess header line for this process.
        Uses the existing header to check whether we require an index
        and (UTC) date appended to the key
        """
        head_key = BASE_VCF_PROCESS_KEY
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
            value=BASE_VCF_PROCESS_LOG.format(
                self.vcfinname, self.run_script, self.arg_str
            ),
        )
        return new_process_line

    def generate_new_increment_header(
        self, existing_line: vcfpy.header.HeaderLine, n: int
    ) -> vcfpy.header.HeaderLine:
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

    def parse_header_add_merge_and_process(
        self, writer_header: vcfpy.Header, max_len: int
    ) -> vcfpy.Header:
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
            for i in range(1, max_len + 1):
                lines_to_add.append(self.generate_new_increment_header(inf_line, i))
        for form_line in format_lines:
            for i in range(1, max_len + 1):
                lines_to_add.append(self.generate_new_increment_header(form_line, i))

        lines_to_add.append(vcfpy.InfoHeaderLine("SPCONF"))
        for new_head_line in lines_to_add:
            writer_header.add_line(new_head_line)

        # Add confidence score info headerline
        ##INFO=<ID=SPCONF,Number=.,Type=String,Description="Smart-Phase confidence scores of the form score_MNV_base_1to2,score_MNV_base_2to3 etc.">
        writer_header.add_info_line(
            vcfpy.OrderedDict(
                [
                    ("ID", "SPCONF"),
                    ("Number", "."),
                    ("Type", "String"),
                    (
                        "Description",
                        "Smart-Phase confidence scores of the form score_MNV_base_1to2,score_MNV_base_2to3 etc.",
                    ),
                ]
            )
        )

        return writer_header

    def merge_snv_to_mnv(
        self, snv_list: List[vcfpy.Record], conf_calls: List
    ) -> vcfpy.Record:
        """
        Merge snvs from list into a single variant and output to VCF
        """
        do_qual = 0
        do_filter = 0
        qual = 0
        filter = []
        if snv_list[0].QUAL:
            do_qual = 1
            LOGGER.warning("Found a QUAL value, output will be a mean of all QUAL.")
        else:
            qual = snv_list[0].QUAL

        do_filter = 0
        filter = []
        for snv in snv_list:
            if snv.FILTER:
                do_filter = 1
                LOGGER.warning(
                    "Found a FILTER value, output will be all FILTERs encountered at all bases."
                )
                break

        chrom = snv_list[0].CHROM
        pos = snv_list[0].POS
        id = []  # list of SNV IDs?
        ref = ""
        alt_str = ""
        info = {"SPCONF": ",".join(map(str, conf_calls))}
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
                # If we want to append qualities to generate a mean
                for f in var.FILTER:
                    filter.append(f)

            if do_qual:
                qual += var.QUAL
            id.append(var.ID[0])

        alt = [vcfpy.Substitution(type_="MNV", value=alt_str)]

        # Make the quality a mean value for MNVs
        if do_qual:
            qual = qual / len(snv_list)

        if "PASS" in filter and all_equal(filter):
            # Check for entries being passes and mark as a single pass
            filter = ["PASS"]

        mnv = vcfpy.Record(
            chrom, pos, id, ref, alt, qual, filter, info, format, calls_dict.values()
        )
        return mnv

    def perform_mnv_merge_to_vcf(self):
        """
        Iterate through VCF records. Outputting a new VCF with
        requested filters removed.
        """
        hom_bed_parsed = None
        if self.bed:
            hom_bed_parsed = parse_homs_bed_to_dict(self.bed)

        (mnvs, max_len) = parse_sphase_output(
            self.spout, self.cutoff, self.exclude_flags, hom_bed_parsed
        )
        reader = self.vcfin
        # Make a copy of the header
        writer_header = reader.header.copy()
        writer_header = self.parse_header_add_merge_and_process(writer_header, max_len)
        writer = vcfpy.Writer.from_path(self.vcfout, writer_header)

        snvs = []
        start_pos_mnv = 0
        start_contig_mnv = ""
        in_mnv = False
        mnv_end_pos = 0
        mnv_conf_calls = None
        for variant in reader:
            if variant.CHROM in mnvs:
                # Start position in an mnv
                if int(variant.POS) in mnvs[variant.CHROM]:
                    in_mnv = True
                    (mnv_end_pos, mnv_conf_calls) = mnvs[variant.CHROM][
                        int(variant.POS)
                    ]
                    start_pos_mnv = int(variant.POS)
                    start_contig_mnv = variant.CHROM
                    snvs.append(variant)
                # In an MNV and waiting for finish
                elif in_mnv and int(variant.POS) <= mnv_end_pos:
                    snvs.append(variant)
                    if int(variant.POS) == mnv_end_pos:
                        mnv_rec = self.merge_snv_to_mnv(snvs, mnv_conf_calls)
                        writer.write_record(mnv_rec)
                        snvs.clear()
                        in_mnv = False
                        mnv_end_pos = 0
                        mnv_conf_calls = None

                else:
                    writer.write_record(variant)
            else:
                writer.write_record(variant)
        writer.close()
