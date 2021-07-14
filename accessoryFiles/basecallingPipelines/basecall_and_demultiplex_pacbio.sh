#!/bin/bash

module load nextflow
module load singularity 

nextflow run -resume -c `pwd`/accessoryFiles/conf/config.nf \
                        `pwd`/pacbio_basecalling_Alleva_et_al_2021.nf \
                        --inPB `pwd`/accessoryFiles/samples/pacbio.details.txt \
                        --pipedir `pwd` \
                        --outdir output 
