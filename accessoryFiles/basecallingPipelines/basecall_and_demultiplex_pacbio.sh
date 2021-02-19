#!/bin/bash
nextflow run -resume -c configPB.nf pacbio_basecalling_Alleva_et_al_2021.nf \
   --inPB    `pwd`/accessoryFiles/samples/pacbio.details.txt \
   --pipedir `pwd`  
