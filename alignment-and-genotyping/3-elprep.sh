#!/bin/bash

library_id=$1
cpus=$2

elprep sfm \
    out/bam/${library_id}.bam \
    out/bam/${library_id}.elprep.bam \
    --nr-of-threads $cpus
    --mark-duplicates \
    --mark-optical-duplicates out/metrics/${library_id}.metrics.txt \
    --sorting-order coordinate \
    --reference indices/bGalGal1b.chr_names.elfasta \
    --haplotypecaller out/vcf/${library_id}.gvcf.gz \
    --intermediate-files-output-type sam \
    --tmp-path $SCRATCH
