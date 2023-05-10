# Genome assembly scripts
This subdirectory contains:
* `assemble.sh`: a SLURM script to assemble the 13 "additional" chicken genomes
  and calculate basic stats
* `assemblies.tsv`: A tsv with two columns: the ID of the individual bird and
  the name of the unaligned bam file containing HiFi reads for that individual.
* `hifiasm.yml`: A conda yml for the environment used to run hifiasm, including
  version numbers.

To actually perform the assemblies, you would need to download the reads, which
are available on SRA. You will likely also need to adapt the SLURM script to
the system you're running it on.
