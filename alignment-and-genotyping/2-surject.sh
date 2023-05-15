#!/bin/bash -e

library_id=$1
cpus=$2

vg surject \
    --threads $cpus \
    --xg-name indices/pangenome.gbz \
    --bam-output \
    out/gam/$library_id.gam \
    | samtools reheader -c "sed s/bGalGal1b#0#//g'" - \
    > out/bam/$library_id.bam
