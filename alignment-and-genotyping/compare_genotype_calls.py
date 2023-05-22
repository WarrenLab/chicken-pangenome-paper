#!/usr/bin/env python3
"""Compare per-sample calls between methods

This is a simple script to compare calls from two methods. As input, it
takes a single vcf(.gz) file that is the output of

`bcftools merge --force-samples caller1_out.vcf caller2_out.vcf`

where caller[1,2]_out.vcf have the same sample names. This script finds
samples that have been given conflict-resolved names by bcftools merge
(e.g., "SAMPLE1" and "2:SAMPLE1") and compares each pair for every
variant, outputting a single number for each variant equal to the number
of sample pairs with the same call divided by the number of sample pairs
where at least one haplotype has a call. If there are no calls, it
outputs -1 to avoid dividing by zero.
"""
import argparse
import gzip
import re
import sys
from typing import IO, Iterable


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "vcf",
        nargs="?",
        help="a vcf file output from `bcftools merge` [STDIN]",
        type=open_gz_or_normal,
        default=sys.stdin,
    )
    return parser.parse_args()


def open_gz_or_normal(vcf_path: str) -> IO:
    """open a file as gzip or normal based on extension

    Given the path to a vcf, figure out if it's gzipped and open it
    properly whether or not it is.
    """
    if vcf_path.endswith(".gz"):
        return gzip.open(vcf_path, "rt")
    else:
        return open(vcf_path, "r")


def parse_vcf(vcf: IO) -> Iterable[float]:
    """parse and analyze a vcf

    Go through a vcf line-by-line, gleaning necessary info from the
    header and then yielding a concordance score for each non-header
    line.
    """
    for line in vcf:
        if line.startswith("##"):
            pass
        elif line.startswith("#"):
            sample_table = make_sample_table(line)
        else:
            yield calc_concordance(sample_table, line)


def make_sample_table(line: str) -> list[tuple[int, int]]:
    """make a sample table from the header line

    Given the last header line of a vcf (i.e., /^#CHROM/), use the
    sample names to make a list of pairs of sample indices where each
    pair is the same sample called by a different method.

    Args:
        line: the text of the last header line of a vcf

    Returns:
        a list of pairs of sample indices
    """
    sample_ids = line.strip().split("\t")[9:]
    table = []
    sample_id_dict: dict[str, int] = {}
    for i, sample_id in enumerate(sample_ids):
        if sample_id.startswith("2:"):
            table.append((sample_id_dict[sample_id[2:]], i))
        else:
            sample_id_dict[sample_id] = i

    return table


gt_seperator_re = re.compile(r"\||\/")


def calc_concordance(sample_table: list[tuple[int, int]], line: str) -> float:
    """calculate the concordance of a single variant

    Calculate the concordance of a single variant line of a vcf. This
    number is equal to the number of sample pairs with the same call
    divided by the number of sample pairs with at least one call.

    Args:
        sample_table: a list of pairs of sample indices indicating that
            the two samples in these columns are the same sample, but
            called by different methods
        line: a plain text line of a vcf

    Returns: a concordance score in the range [0,1], or -1 if there are
        no variant calls on this line
    """
    fields = line.strip().split("\t")
    genotype_field_index = fields[8].split(":").index("GT")

    calls = [s.split(":")[genotype_field_index] for s in fields[9:]]
    num_calls, num_agreements = 0, 0
    for giraffe_call_index, minimap_call_index in sample_table:
        giraffe_call = tuple(gt_seperator_re.split(calls[giraffe_call_index]))
        minimap_call = tuple(gt_seperator_re.split(calls[minimap_call_index]))

        if not (giraffe_call == (".", ".") and minimap_call == (".", ".")):
            num_calls += 1
            if giraffe_call == minimap_call:
                num_agreements += 1

    if num_calls == 0:
        return -1
    else:
        return num_agreements / num_calls


def main():
    """main method of program"""
    args = parse_args()
    for concordance in parse_vcf(args.vcf):
        print(concordance)


if __name__ == "__main__":
    main()
