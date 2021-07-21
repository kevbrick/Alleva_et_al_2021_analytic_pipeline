#!/bin/bash
module load nextflow
module load singularity

if [ ! -f "accessoryFiles/hg38_genome/genome.fa" ]; then 
	echo "Getting hg38 genome FASTA ..."
	bash accessoryFiles/hg38_genome/getGenome.sh
fi

nextflow run Alleva_et_al_2021_analytic_pipeline.nf -c config_Analyses.nf -resume modest_leakey\
  --pipedir      `pwd` \
  --ontfa        `pwd`/fasta/nanopore \
  --pbfa         `pwd`/fasta/pacbio \
  --outdir       output \
  --accessorydir `pwd`/accessoryFiles/ \
  --genomefa     `pwd`/accessoryFiles/hg38_genome/genome.fa \
  --hs           `pwd`/accessoryFiles/DSBhotspots/ \
  --npeaks       1000 \
  --centerSz     1000 \
  --bigwinSize   2000 \
  --RM           true

