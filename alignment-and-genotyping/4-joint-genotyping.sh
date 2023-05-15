#!/bin/bash -e

gatk CombineGVCFs \
    -R indices/bGalGal1b.chr_names.fa \
    $(for vcf in out/vcf/*.gvcf.gz; do
        echo "--variant $vcf"
    done | tr '\n' ' ') \
    -O out/combined.gvcf.gz \
    --tmp-dir $SCRATCH

gatk GenotypeGVCFs \
    -R indices/bGalGal1b.chr_names.fa \
    -V out/combined.gvcf.gz \
    -O out/joint_genotyped.vcf.gz
