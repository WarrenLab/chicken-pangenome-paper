# Minigraph-cactus pangenome creation

We used a nextflow pipeline to make the minigraph-cactus pangenome. To run this
pipeline, we used the command

```bash
nextflow run WarrenLab/minigraph-cactus-nf \
    --seq-file seqFile.tsv \
    --reference bGalGal1b \
    --chromosomes-file chroms.txt
```

The files `seqFile.tsv` and `chroms.txt` are inside this directory. The
pipeline code and instructions for configuring to run on your system are in a
[separate repository](https://github.com/WarrenLab/minigraph-cactus-nf).
