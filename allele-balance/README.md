# Calculating reference bias based on allele balance

In this folder, there are bash scripts and a c program for calculating
reference bias based on allele balance and heterozygous SNV sites. This
analysis is modeled after [Günther et al. (2019)][guenther2019], in which a
site is considered heterozygous if it has coverage of at least 10x and the
minor allele frequency is at least 25%. The C program looks for sites like this
in pileup output, and then outputs the fraction of reads with an alternate
allele at each one.

## Running

First, compile the C program: `gcc calculate_allele_balance.c`. Then, for each
bam, run `1-prepare-bam.sh` (arguments are input bam, output bam, and number of
cpus) to sort and filter the bam in preparation for the analysis, and
`2-calc-allele-balance.sh` (arguments are input bam and outfile) to do the
pileup and calculate the alternate allele frequency at each het site.

Finally, tabulate the output and calculate summary statistics by running
`3-tabulate.sh`.

[guenther2019]: https://doi.org/10.1371/journal.pgen.1008302
