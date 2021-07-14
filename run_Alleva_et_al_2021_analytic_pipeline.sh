#!/bin/bash
module load nextflow
module load singularity

nextflow run Alleva_et_al_2021_analytic_pipeline.nf -c config_Analyses.nf -profile slurm \
  --pipedir      `pwd` \
  --ontfa        `pwd`/fasta/nanopore \
  --pbfa         `pwd`/fasta/pacbio \
  --bonito1dfa   `pwd`/fasta/bonito1d \
  --outdir       output \
  --accessorydir `pwd`/accessoryFiles/ \
  --genomefa     `pwd`/hg38_genome/genome.fa \
  --hs           `pwd`/accessoryFiles/DSBhotspots/ \
  --npeaks       1000 \
  --centerSz     1000 \
  --bigwinSize   2000 \
  --RM           true

