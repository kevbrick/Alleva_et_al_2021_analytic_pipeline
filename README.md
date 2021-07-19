# Alleva et al 2021: Analytic pipeline
Nextflow-based analytic pipeline for long read genotyping of human PRDM9 (Alleva et al., 2021). This pipeline will re-do all analyses and make the figures from Alleva et al. 2021. 

For a stand-alone app to genotype PRDM9 from a long-read FASTA file, see:

`https://github.com/kevbrick/genotype_prdm9_LR`

## Requirements:
* Nextflow 20.10.0+
* Singularity 3.6.4+

### PREREQUISITES: 
The pipeline is built in nextflow (nextflow.io) using singularity. It has been tested using nextflow/20.10.0 and singularity/3.6.4. Earlier versions of nextflow will not work.

### ACCESSORY DATA: 
Post-basecalling FASTA files for all samples are needed to run the pipeline. These can be obtained from:

`https://hpc.nih.gov/~brickkm/alleva2021/alleva2021.tar.gz` 

Extract this tar in the folder that contains the Alleva_et_al_2021_analytic_pipeline.nf script.

### RUN: 
Clone the repo. From the source folder, run:

`bash run_Alleva_et_al_2021_analytic_pipeline.sh`
