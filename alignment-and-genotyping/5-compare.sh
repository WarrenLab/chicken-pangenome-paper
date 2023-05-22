#!/bin/bash

bcftools isec -Oz --write-index \
    -c some \
    -p isec_out \
    giraffe.final.vcf.gz \
    minimap.final.vcf.gz

bcftools merge -Oz \
    -m all \
    --force-samples \
    isec_out/0002.vcf.gz \
    isec_out/0003.vcf.gz \
    > merged.vcf.gz

echo "Giraffe private SNPS: $(bcftools view -i 'QUAL>=10' -H -v snps \
    isec_out/0000.vcf.gz | wc -l)"
echo "Giraffe private indels: $(bcftools view -i 'QUAL>=10' -H -v indels \
    isec_out/0000.vcf.gz | wc -l)"

echo "Minimap private SNPS: $(bcftools view -i 'QUAL>=10' -H -v snps \
    isec_out/0001.vcf.gz | wc -l)"
echo "Minimap private indels: $(bcftools view -i 'QUAL>=10' -H -v indels \
    isec_out/0001.vcf.gz | wc -l)"

echo "Common SNPS: $(bcftools view --threads $SLURM_CPUS_PER_TASK -i 'QUAL>=10' \
    -H -v snps merged.vcf.gz | wc -l)"
echo "Common SNPS: $(bcftools view --threads $SLURM_CPUS_PER_TASK -i 'QUAL>=10' \
    -H -v indels merged.vcf.gz | wc -l)"
