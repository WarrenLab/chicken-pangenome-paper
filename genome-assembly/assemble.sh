#!/bin/bash -e
#SBATCH -a 1-13
#SBATCH --cpus-per-task 14
#SBATCH --mem 150G
#SBATCH -t 2-00:00:00
#SBATCH -o logs/assemble.%a.out
#SBATCH -e logs/assemble.%a.err

eval "$(conda shell.bash hook)"
conda activate hifiasm

assembly_name=$(head -n $SLURM_ARRAY_TASK_ID assemblies.tsv | tail -n 1 | cut -f1)
hifi_bam_paths=$(head -n $SLURM_ARRAY_TASK_ID assemblies.tsv | tail -n 1 | cut -f2)

hifi_fastq_path="out/${assembly_name}/reads.fastq"

mkdir -p out/${assembly_name}

echo "Converting to fastq..."
for bam in $hifi_bam_paths; do
    samtools fastq -@ $SLURM_CPUS_PER_TASK $bam
done > $hifi_fastq_path

awk 'NR%4==2{c+=length($1)} END {print "total bases: " c}' $hifi_fastq_path \
    > ${hifi_fastq_path}.stats

hifiasm $hifi_fastq_path \
    -o out/${assembly_name}/${assembly_name} \
    -t $SLURM_CPUS_PER_TASK

echo "Done assembling..."
rm $hifi_fastq_path

awk '/^S/ { print ">" $2 "\n" $3 }' \
    out/${assembly_name}/${assembly_name}.bp.hap1.p_ctg.gfa \
    > out/${assembly_name}/${assembly_name}.hap1.p_ctg.fa
awk '/^S/ { print ">" $2 "\n" $3 }' \
    out/${assembly_name}/${assembly_name}.bp.hap2.p_ctg.gfa \
    > out/${assembly_name}/${assembly_name}.hap2.p_ctg.fa
awk '/^S/ { print ">" $2 "\n" $3 }' \
    out/${assembly_name}/${assembly_name}.bp.p_ctg.gfa \
    > out/${assembly_name}/${assembly_name}.p_ctg.fa

assembly-stats out/${assembly_name}/${assembly_name}.hap1.p_ctg.fa \
    > out/${assembly_name}/${assembly_name}.hap1.p_ctg.stats
assembly-stats out/${assembly_name}/${assembly_name}.hap2.p_ctg.fa \
    > out/${assembly_name}/${assembly_name}.hap2.p_ctg.stats
assembly-stats out/${assembly_name}/${assembly_name}.p_ctg.fa \
    > out/${assembly_name}/${assembly_name}.p_ctg.stats
