#!/bin/bash -e

cpus=$1

vg paths -F -Q 'bGalGal1b#0#' -x pangenome.gbz > bGalGal1b.paths.fa
minimap2 -ax sr -t $cpus --secondary=no \
    bGalGal1b.paths.fa \
    simulated_1M.interleaved.fastq \
    | awk '/^@/ || $3 != "*"' \
    > mapped.no_unmapped.minimap.sam

vg inject -t $cpus \
    -x pangenome.gbz mapped.minimap.no_unmapped.sam \
    > mapped.minimap.gam

vg annotate -t $cpus \
    -a mapped.minimap.gam \
    -m -x pangenome.gbz \
    | vg gamcompare -t $cpus \
    -s -r 100 -T - simulated_1M.gam \
    > minimap_scores.tsv

echo -e "min_mapQ\treads_aligned\treads_aligned_correctly" \
    > minimap_cumulative_rates.tsv
awk 'NR>1 {aligned[$2]++} NR>1&&$1==1{correct[$2]++}
     END{OFS="\t"; for (i in aligned) print i, aligned[i], correct[i]}' \
    minimap_scores.tsv \
    | sort -nrk1,1 | awk '{OFS="\t"; cum_aligned+=$2; cum_correct+=$3;
         print $1, cum_aligned, cum_correct}' \
    >> minimap_cumulative_rates.tsv
