"""Calculate how much sequence each sample adds to the graph

Given a pangenome graph in GFA format, find how much new sequence each
sample contributes to the graph. To be exact, the algorithm is as
follows:

1. Start with either an empty graph or a graph containing a single
   sample, depending on arguments.
2. Go through all the samples not yet in the graph and determine how
   much new sequence each one would add to the graph.
3. Choose the sample that would add the most new sequence to the graph,
   report it, and add it to the graph.
4. Repeat 2-3 until all sequences are in the graph.

The output is a TSV where the samples are reported in descending order
of how much sequence they add to the graph, with columns sample name and
length of new sequence added.
"""
import argparse
import gzip
import re
from collections import defaultdict
from dataclasses import dataclass
from sys import stderr
from typing import Iterable, Union

WALK_DELIMITER = re.compile(r">|<")


@dataclass
class Segment:
    """Minimal representation of S-line

    For S-lines, we only need to know the node ID and sequence length
    for the purposes of this script, so this is a class that can parse
    the text of an S-line and store those two pieces of information.
    """

    node_id: int
    seq_length: int

    def __init__(self, line: str):
        splits = line.strip().split("\t")
        self.node_id = int(splits[1])
        self.seq_length = len(splits[2])


@dataclass
class Walk:
    """Minimal representation of W-line

    For W-lines, we only need to know what sample it is associated with
    and the IDs of the nodes it traverses, so this class parses the text
    of a W-line and stores these two pieces of information.
    """

    sample_id: str
    nodes: list[int]

    def __init__(self, line):
        # TODO: add a flag to consider the haplotype number as part of the
        # sample ID
        splits = line.strip().split("\t")
        self.sample_id = splits[1]
        self.nodes = list(map(int, WALK_DELIMITER.split(splits[6])[1:]))


@dataclass
class GraphInfo:
    """info about the current state of the graph

    This is a perhaps awkward container for storing two dictionaries
    containing information about the current state of the graph, and
    modifying them as paths iteratively get removed.
    """

    node_lengths: dict[int, int]
    """maps node IDs to length of sequence in the node"""
    sample_node_sets: dict[str, set[int]]
    """maps sample IDs to set of all node IDs used by that sample"""

    def __str__(self):
        return str(self.node_lengths) + "\n" + str(self.sample_node_sets)

    def remove_sample(self, sample_to_remove: str):
        """remove a sample from the graph

        Remove a sample from the sample_node_sets dictionary, and also
        remove all nodes used by this sample from the sets of nodes
        used by all other samples, so that the sample node set for each
        remaining sample contains only the nodes that are unique to it
        compared to all the samples that have been removed so far.
        """
        nodes_to_remove = self.sample_node_sets[sample_to_remove]
        del self.sample_node_sets[sample_to_remove]
        for sample in self.sample_node_sets.keys():
            self.sample_node_sets[sample] -= nodes_to_remove

    def find_longest_sample(self) -> tuple[str, int]:
        """find the sample containing the most unique sequence

        Calculate the amount of unique sequence each sample adds to the
        graph, and then return the name and unique sequence length of
        the sample containing the most unique sequence.
        """
        unique_lengths: list[tuple[str, int]] = []
        for sample, nodes in self.sample_node_sets.items():
            unique_lengths.append((sample, sum(self.node_lengths[n] for n in nodes)))
            print("\t".join(map(str, [sample, nodes, unique_lengths[-1][1]])))
        return max(unique_lengths, key=lambda t: t[1])

    def empty(self) -> bool:
        return len(self.sample_node_sets) == 0


def build_graphinfo(gfa_iter: Iterable[Union[Segment, Walk]]) -> GraphInfo:
    """build a GraphInfo object out of a GFA file"""
    segment_lengths: dict[int, int] = {}
    sample_node_sets: dict[str, set[int]] = defaultdict(set)
    for entry in gfa_iter:
        if isinstance(entry, Segment):
            segment_lengths[entry.node_id] = entry.seq_length
        elif isinstance(entry, Walk):
            sample_node_sets[entry.sample_id].update(entry.nodes)
    return GraphInfo(segment_lengths, sample_node_sets)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-f",
        "--first-haplotype",
        help="remove this haplotype first,"
        " instead of the one with the most unique sequence",
    )
    parser.add_argument("gfa", type=parse_gfa, help="gfa(.gz) file to parse")
    return parser.parse_args()


def parse_gfa(gfa_path: str) -> Iterable[Union[Segment, Walk]]:
    """parse a GFA, yielding segments and walks

    Parse a GFA file. Yields a Segment for S-lines and a Walk for
    W-lines, but nothing for any other lines. Can read plain text or
    gzipped GFA files.
    """
    if gfa_path.endswith(".gz"):
        gfa_file = gzip.open(gfa_path, "rt")
    else:
        gfa_file = open(gfa_path, "r")

    for line in gfa_file:
        if line.startswith("S"):
            yield Segment(line)
        elif line.startswith("W"):
            yield Walk(line)


def main():
    args = parse_args()

    print("Reading graph...", file=stderr)
    graph_info = build_graphinfo(args.gfa)

    if args.first_haplotype:
        print(f"Removing sample {args.first_haplotype} from graph...", file=stderr)
        graph_info.remove_sample(args.first_haplotype)

    i = 0
    while not graph_info.empty():
        print(f"------ROUND {i+1}------", file=stderr)
        most_unique_sample, unique_seq_length = graph_info.find_longest_sample()
        print(f"{most_unique_sample}\t{unique_seq_length}")
        graph_info.remove_sample(most_unique_sample)
        i += 1


if __name__ == "__main__":
    main()
