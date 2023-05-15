#!/bin/bash -e

library_id=$1
cpus=$2

minimap2 -t $cpus -ax sr \
    -R "@RG\tID:$library_id\tSM:$library_id" \
    indices/bGalGal1b.chr_names.idx \
    reads/${library_id}_R1.fastq.gz \
    reads/${library_id}_R2.fastq.gz \
    | samtools view -bh - \
    > out/bam/$library_id.bam
