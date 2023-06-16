# chicken-pangenome-paper

Scripts used to perform analyses in Rice et al. (2023) "A pangenome graph
reference of 30 chicken genomes allows genotyping of large and complex
structural variants."

## Introduction

This repository contains scripts used to perform analyses in the paper "A
pangenome graph reference of 30 chicken genomes allows genotyping of large and
complex structural variants.". It is divided into the following subdirectories,
each with their own README:

* `genome-assembly`: scripts used to assemble genomes to use to create the
  pangenome
* `pangenome-creation`: scripts used to create both the minigraph-cactus and
  PGGB pangenomes
* `pangenome-analysis`: scripts for doing some calculations about the graph and
  its structure
* `alignment-and-genotyping`: scripts used to align reads to the pangenome and
  genotype based on those alignments
* `simulation`: scripts used to create simulated reads, align them to the
  pangenome with different methods, and compare the accuracy of the different
  alignments
* `allele-balance`: looking at reference bias based on allele balance 
