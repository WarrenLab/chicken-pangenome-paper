#!/bin/bash -e
#
library_id=$1
cpus=$2

vg giraffe -t $cpus -p \
    -Z indices/pangenome.gbz \
    -m indices/pangenome.min \
    -d indices/pangenome.dist \
    -f reads/${library_id}_R1.fastq.gz \
    -f reads/${library_id}_R2.fastq.gz \
    -N $library_id -R $library_id \
    > out/gam/${library_id}.gam
