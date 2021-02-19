#!/bin/bash
nextflow run -c configONT.nf nanopore_basecalling_Alleva_et_al_2021.nf \
   --inONT   `pwd`/accessoryFiles/samples/ont.details.txt \
   --pipedir `pwd` \
   --outdir dmx_nanopore
