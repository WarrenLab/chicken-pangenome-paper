#!/usr/bin/env python3
"""
Calculate the allele balance at each heterozygous site in order to
estimate reference bias.
"""

import argparse
import sys
from collections import Counter
from typing import IO


def parse_base_field(base_field: str) -> float:
    """calculate allele balance from base field

    Parse the base field of a pileup. If the position is a
    putative heterozygous site, defined as a site with coverage >= 10
    and minor allele frequency >= 0.25, return the reference allele
    proportion. Otherwise, return -1.

    N.B. this does not work with base fields containing [+-$^<>], so
    mpileup must be run with options
    * --no-output-ends
    * --no-output-ins --no-output-ins (yes, twice)
    * --no-output-del --no-output-del (yes, twice)
    """

    ref_reads_count = 0
    nonref_reads: Counter[str] = Counter()
    for char in base_field:
        if char == "," or char == ".":
            ref_reads_count += 1
        else:
            nonref_reads[char.upper()] += 1

    nonref_base, nonref_base_count = nonref_reads.most_common(1)[0]

    total_reads = ref_reads_count + nonref_base_count
    ref_allele_frequency = ref_reads_count / total_reads
    nonref_allele_frequency = nonref_base_count / total_reads

    if (
        total_reads / len(base_field)
        and min(ref_allele_frequency, nonref_allele_frequency) > 0.25
    ):
        return nonref_allele_frequency

    return -1


def parse_pileup(infile: IO):
    """Parse a pileup file and calculate reference frequencies

    Parse a pileup line by line. For every putative heterozygous site,
    defined here as a site with coverage of at least 10 where more than
    25% of the reads contain the minor allele, yield the reference
    allele frequency.
    """
    for line in infile:
        fields = line.strip().split("\t")
        depth = int(fields[3])
        if depth >= 10 and depth <= 100:
            het_site_ref_allele_frequency = parse_base_field(fields[4])
            if het_site_ref_allele_frequency != -1:
                yield het_site_ref_allele_frequency


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "pileup",
        type=lambda p: parse_pileup(open(p, "r")),
        default=parse_pileup(sys.stdin),
        help="samtools pileup output to read [STDIN]",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    for het_site_ref_allele_frequency in args.pileup:
        print(het_site_ref_allele_frequency)


if __name__ == "__main__":
    main()
