## Pipelines for Nanopore and Pacbio basecalling for Alleva et al 2021 (Genotyping PRDM9 in humans)

### PREREQUISITES: 
These pipelines are built in nextflow (nextflow.io) and use either singularity or coda for package management. It has been tested using nextflow/20.10.0 and singularity/3.6.4. Earlier versions of nextflow will not work.

### ACCESSORY DATA: 
The raw Pacbio CCS BAM files and Nanopore FAST5 files are available upon request. After discussion with the SRA / GEO, these files have not been included in public repositorys. For more information, please contact the authors. 

### RUN FOR NANOPORE: 
Clone the repo. Modify the accessoryFiles/samples/ont.details.txt file to specify the folders for the FAST5 files for each ONT run. Then, from the source folder, run:

`bash basecall_and_demultiplex_nanopore.sh`

### RUN FOR PACBIO: 
Clone the repo. Modify the accessoryFiles/samples/pb.details.txt file to specify the folders for the Pacbio-CCS BAM file for each run. Then, from the source folder, run:

`bash basecall_and_demultiplex_pacbio.sh`
