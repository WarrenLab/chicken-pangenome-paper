#!/bin/bash -e

vg sim --progress \
    --xg-name pangenome.gbz \
    --num-reads 1000000 \
    --read-length 150 \
    --align-out \
    --random-seed 12345 \
    --sub-rate 0.0024 \
    --indel-rate 0.00029 \
    --frag-len 570 \
    --frag-std-dev 165 \
    > simulated_1M.gam

vg view -X -a simulated_1M.gam > simulated_1M.interleaved.fastq
