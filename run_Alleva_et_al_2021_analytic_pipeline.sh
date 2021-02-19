#!/bin/bash
nextflow run -c config_Analyses.nf Alleva_et_al_2021_analytic_pipeline.nf \
  --pipedir      `pwd` \
  --ontfa        `pwd`/fasta/nanopore \
  --pbfa         `pwd`/fasta/pacbio \
  --outdir       `pwd`/output \
  --accessorydir `pwd`/accessoryFiles/ \
  --genomefa     `pwd`/hg38_genome/genome.fa \
  --hs           `pwd`/DSBhotspots/ \
  --npeaks       1000 \
  --centerSz     1000 \
  --bigwinSize   2000 \
  --RM           true

