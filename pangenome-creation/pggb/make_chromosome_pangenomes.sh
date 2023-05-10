#!/bin/bash -e
#SBATCH -a 1-41
#SBATCH --mem 150G
#SBATCH --cpus-per-task 13
#SBATCH -t 2-00:00:00
#SBATCH -o logs/pggb.%a.out
#SBATCH -e logs/pggb.%a.err

eval "$(conda shell.bash hook)"
conda activate pggb

line_number=$((SLURM_ARRAY_TASK_ID + 1))
bGalGal5_chr=$(head -n $line_number final_params.csv | tail -n 1 | cut -f1 -d,)
bGalGal1b_chr=$(head -n $line_number final_params.csv | tail -n 1 | cut -f2 -d,)
bGalGal1w_chr=$(head -n $line_number final_params.csv | tail -n 1 | cut -f3 -d,)
bGalGal4_chr=$(head -n $line_number final_params.csv | tail -n 1 | cut -f4 -d,)
huxuT2T_chr=$(head -n $line_number final_params.csv | tail -n 1 | cut -f5 -d,)
segment_length=$(head -n $line_number final_params.csv | tail -n 1 | cut -f6 -d,)
map_pct_id=$(head -n $line_number final_params.csv | tail -n 1 | cut -f7 -d,)
min_match_len=$(head -n $line_number final_params.csv | tail -n 1 | cut -f8 -d,)

chromosome_name=$bGalGal5_chr
mkdir -p $chromosome_name
cd $chromosome_name

refname="bGalGal1b"; refChr=$bGalGal1b_chr
samtools faidx ../../references/${refname}.fa ${refChr} | \
    sed "s/^>${refChr}/>${refname}#${refChr}/" > ${chromosome_name}.fa

# bGalGalw does not have sex chromosomes
refname="bGalGal1w"; refChr=$bGalGal1w_chr
if [ ! -z "$refChr" ]; then
    samtools faidx ../../references/${refname}.fa ${refChr} | \
        sed "s/^>${refChr}/>${refname}#${refChr}/" >> ${chromosome_name}.fa
fi

refname="bGalGal4"; refChr=$bGalGal4_chr
samtools faidx ../../references/${refname}.fa ${refChr} | \
    sed "s/^>${refChr}/>${refname}#${refChr}/" >> ${chromosome_name}.fa
refname="bGalGal5"; refChr=$bGalGal5_chr
samtools faidx ../../references/${refname}.fa ${refChr} | \
    sed "s/^>${refChr}/>${refname}#${refChr}/" >> ${chromosome_name}.fa
refname="huxuT2T"; refChr=$huxuT2T_chr
samtools faidx ../../references/${refname}.fa ${refChr} | \
    sed "s/^>${refChr}/>${refname}#${refChr}/" >> ${chromosome_name}.fa

bgzip ${chromosome_name}.fa
samtools faidx ${chromosome_name}.fa.gz

# run pggb
pggb \
    -i ${chromosome_name}.fa.gz \
    -n $(wc -l ${chromosome_name}.fa.gz.fai) \
    -t $SLURM_CPUS_PER_TASK \
    -s $segment_length \
    -p $map_pct_id \
    -k $min_match_len \
    -G 3079,3559 \
    -V 'bGalGal5:#'

# make a better-looking 2D visualization
odgi draw -i *.smooth.final.og -c *.smooth.final.og.lay -C -w1000 \
    -p ${chromosome_name}.2d.png

# make a better-looking 1D visualization
odgi sort -H <(echo "bGalGal5#${chromosome_name}") -Y \
    -i *.smooth.final.og -o ${chromosome_name}.smooth.final.sorted.og \
    -t $SLURM_CPUS_PER_TASK
odgi viz -t $SLURM_CPUS_PER_TASK \
    -i ${chromosome_name}.smooth.final.sorted.og \
    -o ${chromosome_name}.1d.png
