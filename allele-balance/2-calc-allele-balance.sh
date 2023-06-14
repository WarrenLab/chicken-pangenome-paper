#!/bin/bash -e

inbam=$1
outfile=$2

samtools mpileup \
    -q 10 \
    --no-output-ends \
    --no-output-ins --no-output-ins \
    --no-output-del --no-output-del \
    --no-BAQ \
    -d 100 \
    -f bGalGal1b.chr_names.fa \
    $infile \
    | ./a.out \
    > $outfile
