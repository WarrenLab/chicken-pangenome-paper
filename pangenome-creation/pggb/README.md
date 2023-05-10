# PGGB pangenome creation

To create the PGGB pangenome, we first created a pangenome for each chromosome
using `make_chromosome_pangenomes.sh`, and then merged them with
`merge_chromosomes.sh`. These are both SLURM scripts, so you will likely need
to modify them to run on your system. The parameters file `final_params.csv` is
also in this directory.

There is also a conda environment export file `pggb.yml` listing the software
used to run this pipeline, including version numbers.
