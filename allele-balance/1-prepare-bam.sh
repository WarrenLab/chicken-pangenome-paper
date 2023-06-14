#!/bin/bash -e

inbam=$1
outbam=$2
cpus=$3

samtools view -bh -q10 -F2304 -@ $cpus \
    $inbam \
    | samtools sort -@ $cpus \
    -T $TMPDIR/$(basename $inbam) - \
    > $outbam
