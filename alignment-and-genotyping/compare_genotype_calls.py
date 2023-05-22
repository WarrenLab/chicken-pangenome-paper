#!/usr/bin/env python3

import argparse
import gzip
import re
from typing import Iterable


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "vcf", help="a vcf file output from `bcftools merge`", type=parse_vcf
    )
    return parser.parse_args()


def parse_vcf(vcf_path: str) -> Iterable[float]:
    for line in gzip.open(vcf_path, "rt"):
        if line.startswith("##"):
            pass
        elif line.startswith("#"):
            sample_table = make_sample_table(line)
        else:
            yield calc_concordance(sample_table, line)


def make_sample_table(line) -> list[tuple[int, int]]:
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
    args = parse_args()
    for concordance in args.vcf:
        print(concordance)


if __name__ == "__main__":
    main()
