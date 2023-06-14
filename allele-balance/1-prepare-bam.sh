#!/bin/bash -e

library_id=$1
aligner=$2
cpus=$3
scratch=/scratch/allele_balance

samtools view -bh -q10 -F2304 -@ $cpus \
    $scratch/bams/${library_id}.${aligner}.bam \
    | samtools sort -@ $cpus \
    -T $scratch/bams/${library_id}.${aligner} - \
    > $scratch/bams/${library_id}.${aligner}.s.bam
