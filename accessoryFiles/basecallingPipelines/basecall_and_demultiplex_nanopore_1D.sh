#!/bin/bash

module load nextflow/20.07.1
module load singularity 

nextflow run -resume -c configONT.nf  nanopore_basecalling_bonito1D_Alleva_et_al_2021.nf \
                        --inONT `pwd`/accessoryFiles/samples/ont.details.txt  \
                        --pipedir `pwd`  \
                        --outdir dmx_nanopore \
                        --fast5groupsize 50

