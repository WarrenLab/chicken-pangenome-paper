#!/bin/bash
#SBATCH --cpus-per-task 14
#SBATCH --mem 100G
#SBATCH -t 1-00:00:00
#SBATCH -o combine_chromosomes.out
#SBATCH -e combine_chromosomes.err

eval "$(conda shell.bash hook)"
conda activate pggb

ls chr*/*.smooth.final.sorted.og > ogs_to_squeeze.txt

odgi squeeze \
    -t $SLURM_CPUS_PER_TASK \
    -f ogs_to_squeeze.txt \
    -o whole_chicken_genome.og

ls --no-color chr*/*.vcf | grep -v chrW | grep -v chrZ > vcfs_to_concat.txt
bcftools concat -f vcfs_to_concat.txt > whole_chicken_genome.vcf
