#!/bin/bash -e

scratch=$SCRATCH/allele_balance

echo -e "library\taligner\tmean\tstderr" > alt_freq_stats.tsv
while read library_id; do
    for aligner in giraffe minimap; do
        mean=$(awk '{x+=$6}END{print x/NR}' \
            $scratch/out/${library_id}.${aligner}.alt_freqs)
        stdev=$(awk -v mean=$mean '{s+=($6-mean)^2} END {print sqrt(s/NR)/sqrt(NR)}' \
            $scratch/out/${library_id}.${aligner}.alt_freqs)
        echo -e "$library_id\t$aligner\t$mean\t$stdev"
    done
done < library_ids.txt >> alt_freq_stats.tsv
