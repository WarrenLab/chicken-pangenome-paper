#!/bin/bash

cpus=$1

vg giraffe -t $cpus -p \
    -Z pangenome.gbz \
    -m pangenome.min \
    -d pangenome.dist \
    -i -f simulated_1M.interleaved.fastq \
    > mapped.giraffe.gam

vg annotate -t $cpus -a mapped.giraffe.gam -m -x pangenome.gbz \
    | vg gamcompare -T -t $cpus -s -r 100 - simulated_1M.gam \
    > giraffe_scores.tsv

echo -e "min_mapQ\treads_aligned\treads_aligned_correctly" \
    > giraffe_cumulative_rates.tsv
awk 'NR>1 {aligned[$2]++} NR>1&&$1==1{correct[$2]++}
     END{OFS="\t"; for (i in aligned) print i, aligned[i], correct[i]}' \
    giraffe_scores.tsv \
    | sort -nrk1,1 | awk '{OFS="\t"; cum_aligned+=$2; cum_correct+=$3;
         print $1, cum_aligned, cum_correct}' \
    >> giraffe_cumulative_rates.tsv
