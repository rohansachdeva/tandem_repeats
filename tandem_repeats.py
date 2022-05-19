#!/usr/bin/env python3
import argparse
import functools
import math
import operator
import os
import shutil
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from subprocess import PIPE, Popen

import pandas as pd
import smart_open
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.SeqUtils import GC

__author__ = "Rohan Sachdeva"
__copyright__ = "Copyright 2021"
__maintainer__ = "Rohan Sachdeva"
__email__ = "rohansach@berkeley.edu"
__status__ = "Development"


def create_output(out_dir):
    try:
        os.makedirs(out_dir)

    except FileExistsError:
        if args.force:
            shutil.rmtree(out_dir)
            os.makedirs(out_dir)

        else:
            print("Out folder exists (pass --force to overwrite)", file=sys.stderr)
            sys.exit(1)


def generate_fa(header, seq):
    fa = ">" + header + "\n" + seq

    return fa


def generate_sequence_groups(in_seq_path, min_sequence_len, min_gc, max_gc):
    for header, seq in SimpleFastaParser(smart_open.open(in_seq_path)):

        if len(seq) >= min_sequence_len:
            gc = GC(seq)

            if gc >= min_gc and gc <= max_gc:
                id_ = header.split()[0]
                seq = seq.upper()

                id_seq_group = (id_, seq)

                yield id_seq_group


def get_shannon_entropy_seq(seq):
    seq = seq.upper()

    nt_to_count = defaultdict(int)
    seq_len = len(seq)

    shannon_entropy_seq = 0

    for nt in seq:
        nt_to_count[nt] += 1

    for nt, nt_count in nt_to_count.items():
        nt_frac = nt_count / seq_len

        nt_entropy = nt_frac * math.log(nt_frac, 2)

        shannon_entropy_seq += nt_entropy

    shannon_entropy_seq = shannon_entropy_seq * -1

    return shannon_entropy_seq


def run_tandem_repeat_match(id_seq, repeat_match_path, min_repeat_len):

    cmd = [repeat_match_path, "-t", "-n", str(min_repeat_len), "/dev/stdin"]

    id_, seq = id_seq

    fa = generate_fa(id_, seq)

    cmd = " ".join(cmd)

    proc = Popen(cmd, stdout=PIPE, stdin=PIPE, stderr=PIPE, shell=True)

    proc_stdout = proc.communicate(input=fa.encode())[0].decode().strip().split("\n")

    return proc_stdout


def get_tandem_repeats(
    id_seq,
    min_repeat_len,
    min_repeat_count,
    min_uniq_nt,
    min_region_len,
    up_len,
    down_len,
    tandem_repeat_match_result,
):
    tandem_repeat_match_result = tandem_repeat_match_result[2:]

    id_, seq = id_seq

    tandem_infos = []
    tandem_repeat_match_result_vals = []
    parent_up_starts_uniq = set()

    for line in tandem_repeat_match_result:
        line = line.strip().split()
        line = [int(col) for col in line]

        tandem_repeat_match_result_vals.append(line)

    try:
        tandem_repeat_match_result_vals = sorted(
            tandem_repeat_match_result_vals, key=operator.itemgetter(1, 2)
        )
    except IndexError:
        print(
            "mummer failed",
            id_,
            tandem_repeat_match_result,
            tandem_repeat_match_result_vals,
            file=sys.stderr,
        )

    repeat_num = 0

    for line in tandem_repeat_match_result_vals:
        parent_up_start, parent_down_start, parent_repeat_len = line

        parent_up_end = parent_up_start + parent_repeat_len
        parent_down_end = parent_down_start + parent_repeat_len

        if (
            parent_up_start < parent_down_start
            and parent_up_start + parent_repeat_len > parent_down_start
            and parent_repeat_len >= min_repeat_len
            and parent_up_start not in parent_up_starts_uniq
        ):

            parent_up_starts_uniq.add(parent_up_start)

            tandem_repeat_len = parent_down_start - parent_up_start
            parent_repeat_len = tandem_repeat_len + parent_repeat_len

            tandem_repeat_count = parent_repeat_len / tandem_repeat_len
            tandem_repeat_count = math.trunc(tandem_repeat_count)

            tandem_region_len = tandem_repeat_len * tandem_repeat_count
            tandem_region_end = parent_up_start + tandem_region_len - 1

            tandem_region_seq = seq[parent_up_start - 1 : tandem_region_end - 1 + 1]
            tandem_repeat_seq = seq[
                parent_up_start - 1 : parent_up_start + tandem_repeat_len - 1
            ]

            tandem_repeat_uniq_nts = set(tandem_repeat_seq)
            tandem_repeat_uniq_nts_count = len(set(tandem_repeat_uniq_nts))

            if (
                tandem_repeat_len >= min_repeat_len
                and tandem_repeat_count >= min_repeat_count
                and tandem_region_len >= min_region_len
                and tandem_repeat_uniq_nts_count >= min_uniq_nt
            ):
                repeat_num += 1

                tandem_repeat_seq_gc = GC(tandem_repeat_seq)
                tandem_repeat_shannon_entropy = get_shannon_entropy_seq(
                    tandem_repeat_seq
                )

                up_region_seq = seq[parent_up_start - 1 - up_len : parent_up_start - 1]
                down_region_seq = seq[tandem_region_end : tandem_region_end + down_len]

                tandem_info = {
                    "repeat_num": repeat_num,
                    "region_start": parent_up_start,
                    "region_end": tandem_region_end,
                    "repeat_length": tandem_repeat_len,
                    "repeat_count": tandem_repeat_count,
                    "region_length": tandem_region_len,
                    "repeat_uniq_nts_count": tandem_repeat_uniq_nts_count,
                    "repeat_uniq_nts": "".join(sorted(tandem_repeat_uniq_nts)),
                    "repeat_gc": round(tandem_repeat_seq_gc, 2),
                    "repeat_shannon_entropy": round(tandem_repeat_shannon_entropy, 2),
                    "repeat_seq": tandem_repeat_seq,
                    "up_region_seq": up_region_seq,
                    "down_region_seq": down_region_seq,
                    "region_seq": tandem_region_seq,
                }

                tandem_infos.append(tandem_info)

    if len(tandem_infos) > 0:
        tandem_result = (id_, seq, tandem_infos)

        return tandem_result


def run_get_tandem_repeats(
    repeat_match_path,
    min_repeat_len,
    min_repeat_count,
    min_uniq_nt,
    min_region_len,
    up_len,
    down_len,
    id_seq,
):
    tandem_repeat_match_result = run_tandem_repeat_match(
        id_seq, repeat_match_path, min_repeat_len
    )

    tandem_repeats_result = get_tandem_repeats(
        id_seq,
        min_repeat_len,
        min_repeat_count,
        min_uniq_nt,
        min_region_len,
        up_len,
        down_len,
        tandem_repeat_match_result,
    )

    return tandem_repeats_result


def parallel_run_get_tandem_repeats(
    procs,
    id_seq_groups,
    repeat_match_path,
    min_repeat_len,
    min_repeat_count,
    min_uniq_nt,
    min_region_len,
    up_len,
    down_len,
):

    partial_run_get_tandem_repeats = functools.partial(
        run_get_tandem_repeats,
        repeat_match_path,
        min_repeat_len,
        min_repeat_count,
        min_uniq_nt,
        min_region_len,
        up_len,
        down_len,
    )

    with ProcessPoolExecutor(procs) as executor:
        execute_result = executor.map(partial_run_get_tandem_repeats, id_seq_groups)

        return execute_result


def write_tandem_output(parallel_tandem_results):

    tandem_info_type_to_vals = defaultdict(list)
    id_to_seq = {}

    out_dir_path = args.out_dir

    out_dir_name = os.path.basename(out_dir_path)

    out_tables_dir_path = os.path.join(out_dir_path, "tables")
    out_seqs_dir_path = os.path.join(out_dir_path, "sequences")
    out_annotations_dir_path = os.path.join(out_dir_path, "annotations")

    out_table_name = ".".join((out_dir_name, "info-per-repeat", "tsv"))
    out_table_path = os.path.join(out_tables_dir_path, out_table_name)

    out_seqs_table_name = ".".join((out_dir_name, "info-per-sequence", "tsv"))
    out_seqs_table_path = os.path.join(out_tables_dir_path, out_seqs_table_name)

    out_bed_name = ".".join((out_dir_name, "bed"))
    out_bed_path = os.path.join(out_annotations_dir_path, out_bed_name)

    out_seqs_name = ".".join((out_dir_name, "sequences-with-repeats", "fa"))
    out_seqs_path = os.path.join(out_seqs_dir_path, out_seqs_name)

    create_output(out_dir_path)

    for dir_path in (out_tables_dir_path, out_seqs_dir_path, out_annotations_dir_path):
        os.mkdir(dir_path)

    with open(out_seqs_path, "w") as out_seqs_path_fh:
        for tandem_result in parallel_tandem_results:
            if tandem_result is not None:
                id_, seq, tandem_infos = tandem_result

                id_to_seq[id_] = seq

                fa = generate_fa(id_, seq)
                out_seqs_path_fh.write(fa + "\n")

                for tandem_info_type_to_val in tandem_infos:
                    tandem_info_type_to_vals["sequence_id"].append(id_)

                    for info_type, val in tandem_info_type_to_val.items():
                        tandem_info_type_to_vals[info_type].append(val)

    df_tandem_info = pd.DataFrame(tandem_info_type_to_vals)
    df_tandem_info.to_csv(out_table_path, sep="\t", index=False)

    with open(out_bed_path, "w") as out_bed_path_fh:
        for (
            sequence_id,
            repeat_num,
            region_start,
            region_end,
            region_length,
            repeat_count,
            repeat_length,
        ) in zip(
            df_tandem_info["sequence_id"],
            df_tandem_info["repeat_num"],
            df_tandem_info["region_start"],
            df_tandem_info["region_end"],
            df_tandem_info["region_length"],
            df_tandem_info["repeat_count"],
            df_tandem_info["repeat_length"],
        ):

            repeat_name = "_".join(("tandem", str(repeat_num)))
            region_start_bed = region_start - 1

            repeat_lengths = ",".join([str(repeat_length)] * repeat_count)

            block_starts = []

            for block_start in range(0, region_length, repeat_length):
                block_starts.append(str(block_start))

            block_starts = ",".join(block_starts)

            bed_line = (
                sequence_id,
                region_start_bed,
                region_end,
                repeat_name,
                0,
                "+",
                region_start_bed,
                region_end,
                "255,192,203",
                repeat_count,
                repeat_lengths,
                block_starts,
            )

            bed_line = "\t".join((str(i) for i in bed_line))

            out_bed_path_fh.write(bed_line + "\n")

    summary_info_to_vals = defaultdict(list)

    for id_, group in df_tandem_info.groupby("sequence_id"):
        seq = id_to_seq[id_]
        seq_len = len(seq)
        seq_gc = GC(seq)

        total_repeat_count = group["repeat_count"].sum()
        total_repeat_count_per_bp = total_repeat_count / seq_len
        total_repeat_count_per_len = total_repeat_count_per_bp / args.bp_norm_len

        total_region_len = group["region_length"].sum()
        total_region_percent = total_region_len / seq_len * 100

        total_region_count = len(group["repeat_count"])

        total_region_count_per_bp = total_region_count / seq_len
        total_region_count_per_len = total_region_count_per_bp / args.bp_norm_len

        summary_info_to_vals["sequence_id"].append(id_)
        summary_info_to_vals["sequence_length"].append(seq_len)
        summary_info_to_vals["sequence_gc"].append(round(seq_gc, 2))

        summary_info_to_vals["total_region_count"].append(total_region_count)
        summary_info_to_vals["total_region_count_per_bp"].append(
            total_region_count_per_bp
        )
        summary_info_to_vals[
            "_".join(("total_region_count_per", str(args.bp_norm_len)))
        ].append(total_region_count_per_len)
        summary_info_to_vals["total_region_length"].append(total_region_len)
        summary_info_to_vals["total_region_percent"].append(total_region_percent)

        summary_info_to_vals["total_repeat_count"].append(total_repeat_count),
        summary_info_to_vals["total_repeat_count_per_bp"].append(
            total_repeat_count_per_bp
        ),
        summary_info_to_vals[
            "_".join(("total_repeat_count_per", str(args.bp_norm_len)))
        ].append(total_repeat_count_per_len)

    df_sequence_info = pd.DataFrame(summary_info_to_vals)

    df_sequence_info.to_csv(out_seqs_table_path, sep="\t", index=False)


def main():
    id_seq_groups = generate_sequence_groups(
        args.in_seqs, args.min_seq_len, args.min_gc, args.max_gc
    )

    tandem_results = parallel_run_get_tandem_repeats(
        args.procs,
        id_seq_groups,
        "repeat-match",
        args.min_repeat_len,
        args.min_repeat_count,
        args.min_uniq_nt,
        args.min_region_len,
        args.up_len,
        args.down_len,
    )

    write_tandem_output(tandem_results)


if __name__ == "__main__":
    argparser = argparse.ArgumentParser()

    argparser.add_argument(
        "-i",
        "--in",
        action="store",
        dest="in_seqs",
        required=True,
        help="Input sequences in FASTA format",
    )

    argparser.add_argument(
        "-l",
        "--min_repeat_len",
        action="store",
        dest="min_repeat_len",
        required=False,
        type=int,
        default=9,
        help="Minimum repeat length to find",
    )

    argparser.add_argument(
        "--min_repeat_count",
        action="store",
        dest="min_repeat_count",
        required=False,
        type=int,
        default=2,
        help="Minimum repeat count to find",
    )

    argparser.add_argument(
        "-m",
        "--min_seq_len",
        action="store",
        dest="min_seq_len",
        required=False,
        type=int,
        default=0,
        help="Minimum length of sequence to look for repeats",
    )

    argparser.add_argument(
        "--min_gc",
        action="store",
        dest="min_gc",
        required=False,
        type=float,
        default=0,
        help="Minimum percent GC to look for repeats",
    )

    argparser.add_argument(
        "--max_gc",
        action="store",
        dest="max_gc",
        required=False,
        type=float,
        default=200,
        help="Maximum percent GC to look for repeats",
    )

    argparser.add_argument(
        "--min_uniq_nt",
        action="store",
        dest="min_uniq_nt",
        required=False,
        type=float,
        default=1,
        help="Minimum unique NTs in each tandem repeat",
    )

    argparser.add_argument(
        "--min_region_len",
        action="store",
        dest="min_region_len",
        required=False,
        type=float,
        default=0,
        help="Minimum length of the region made of tandems",
    )

    argparser.add_argument(
        "--up_len",
        action="store",
        dest="up_len",
        required=False,
        type=int,
        default=50,
        help="Minimum length of the region upstream of the repeat region",
    )

    argparser.add_argument(
        "--down_len",
        action="store",
        dest="down_len",
        required=False,
        type=int,
        default=50,
        help="Minimum length of the region downstream of the repeat region",
    )

    argparser.add_argument(
        "--bp_norm_len",
        action="store",
        dest="bp_norm_len",
        required=False,
        type=int,
        default=100000,
        help="BP to normalize to in summaries",
    )

    argparser.add_argument(
        "-p",
        "--procs",
        action="store",
        dest="procs",
        required=False,
        type=int,
        default=1,
        help="Processes to use",
    )

    argparser.add_argument(
        "--force",
        action="store_true",
        dest="force",
        required=False,
        help="Number of split FASTAs to create",
    )

    argparser.add_argument(
        "-o",
        "--out",
        action="store",
        dest="out_dir",
        required=True,
        help="Output directory path",
    )

    args = argparser.parse_args()

    main()
