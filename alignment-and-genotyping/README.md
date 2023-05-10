# Alignment and genotyping

This directory contains scripts used for aligning short reads to the pangenome
and genotyping the individuals according to these alignments.

The scripts are numbered in order: `1-giraffe.sh`, `2-surject.sh`,
`3-genotype.sh`. Each must be run on each library, with two arguments: library
ID and number of threads. For example, to run the whole pipeline on sample
#47 with 10 threads:

```bash
./1-giraffe.sh 47_Standard_Oriental_Phoenix_Silver 10
./2-surject.sh 47_Standard_Oriental_Phoenix_Silver 10
./3-genotype.sh 47_Standard_Oriental_Phoenix_Silver 10
```

You'll need the fastqs in `reads/` and the minigraph-cactus pangenome indices
(.gbz, .min, and .dist) in `indices/`. It will output graph alignments for the
library to `out/bam/${library_id}.gam`, a bam of the alignments to `out/gam/`,
and a gvcf to `out/vcf/`.

We used a SLURM wrapper script to run all the samples though each step
in parallel, but the best way to do that will depend on your system.

To do linear alignment instead of pangenome alignment (for comparison),
substitute `1-minimap.sh` for `1-giraffe.sh`.

## Versions
We used the following software versions in this pipeline:
* vg v1.47.0 (binary from github release)
* samtools v1.17 with htslib 1.17 (compiled from github source release)
* elprep v5.1.2 (binary from github release)
* minimap2 v2.26 (binary from github release)
