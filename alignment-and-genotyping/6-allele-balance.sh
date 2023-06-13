#!/bin/bash -e

library_id=$1

samtools mpileup \
    -q 10 \
    --no-output-ends \
    --no-output-ins --no-output-ins \
    --no-output-del --no-output-del \
    --no-BAQ \
    -d 100 \
    -f bGalGal1b.chr_names.fa \
    out/bam/${library_id}.bam \
    | python calculate_allele_balance.py \
    > out/allele_balance/${library_id}.het_site_ref_freqs
