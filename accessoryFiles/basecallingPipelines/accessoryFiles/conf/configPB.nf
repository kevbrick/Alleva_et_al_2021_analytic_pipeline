params.pipedir="."
//params.outdir="${launchDir}/${params.outdir}"
params.accessorydir = "${params.pipedir}/accessoryFiles/"

timeline {
  enabled = false
  file = "${params.outdir}/nextflow/prdm9_pacbio_1kg_timeline.html"
}

report {
  enabled = false
  file = "${params.outdir}/nextflow/prdm9_pacbio_1kg_report.html"
}

trace {
  enabled = false
  file = "${params.outdir}/nextflow/prdm9_pacbio_1kg_trace.txt"
}

manifest {
  description = 'Human 1KG PRDM9 genotyping pipeline. Author: Kevin Brick.'
}

singularity.enabled = true
singularity.autoMounts = true
singularity.envWhitelist='https_proxy,http_proxy,ftp_proxy,DISPLAY,NXF_GENOMES'
singularity.runOptions=" -B \$NXF_GENOMES -B ${params.accessorydir} -B ${params.pipedir} -B ${launchDir}/${params.outdir} --nv"
executor.$slurm.pollInterval = '1 min'
executor.$slurm.queueStatInterval = '5 min'

process{
  withName:splitFAST5s{
    cpus = { 2 }
    memory = { 4.GB }
    time = { 0.5.h * task.attempt }
  }

  withName:getReadPairs_ONT{
    time           = { 12.h * task.attempt }
    //container      = "docker://kevbrick/ont_guppybonito:1.0"
    module = 'guppy/3.1.5'
    clusterOptions = ' --partition=gpu --cpus-per-task=6 --mem=16g --gres=lscratch:500,gpu:p100:1 '
  }

  withName:splitIndividualsForBonito{
    time           = { 3.h * task.attempt }
    container      = "docker://kevbrick/ont_guppybonito:1.0"
    clusterOptions = ' --partition=gpu --cpus-per-task=6 --mem=16g --gres=lscratch:500,gpu:p100:1 '
  }

  withName:splitFAST5sByIndividual{
    time           = { 8.h * task.attempt }
    container      = "docker://kevbrick/ont_guppybonito:1.0"
    clusterOptions = ' --partition=gpu --cpus-per-task=6 --mem=16g --gres=lscratch:500,gpu:p100:1 '
  }

  withName:bonitoBasecall_ONT{
    time           = { 8.h * task.attempt }
    container      = "docker://kevbrick/ont_guppybonito:1.0"
    clusterOptions = ' --partition=gpu --cpus-per-task=6 --mem=16g --gres=lscratch:500,gpu:p100:1 '
    //conda          = "${params.accessorydir}/conda/bonito_environment.yml"
  }

  withName:split_ONT1d_by_barcode{
    cpus = { 2 }
    memory = { 8.GB }
    time = { 0.5.h * task.attempt }
  }

  withName:getNamedExptChannelPB{
    cpus = { 2 }
    memory = { 8.GB }
    time = { 0.5.h * task.attempt }
    container     = "docker://kevbrick/pacbioccstools:1.0"
    //process.conda = "${params.accessorydir}/conda/environment.yml"
  }

  withName:get_known_human_prdm9_data{
    cpus = { 2 }
    memory = { 8.GB }
    time = { 0.5.h * task.attempt }
  }

  withName:split_pacbio_bam{
    cpus = { 2 }
    memory = { 8.GB }
    time = { 4.h * task.attempt }
    container  	  = "docker://kevbrick/pacbioccstools:1.0"
    //process.conda = "${params.accessorydir}/conda/environment.yml"
  }

  withName:call_ccs{
    cpus = { 16 }
    memory = { 32.GB }
    time = { 16.h * task.attempt }
    container  	  = "docker://kevbrick/pacbioccstools:1.0"
    //process.conda = "${params.accessorydir}/conda/environment.yml"
  }

  withName:demultiplex{
    cpus = { 16 }
    memory = { 32.GB }
    time = { 16.h * task.attempt}
    container  	  = "docker://kevbrick/pacbioccstools:1.0"
    //process.conda = "${params.accessorydir}/conda/environment.yml"
  }

  withName:gatherCCSDemuxFAs{
    cpus = { 2 }
    memory = { 8.GB }
    time = { 3.h * task.attempt}
    container  	  = "docker://kevbrick/pacbioccstools:1.0"
    //process.conda = "${params.accessorydir}/conda/environment.yml"
  }

  withName:merge_raw_fa{
    cpus = { 1 }
    memory = { 8.GB }
    time = { 2.h * task.attempt}
  }

  withName:dmxBAMtoFA{
    cpus = { 2 }
    memory = { 16.GB }
    time = { 4.h * task.attempt}
    container  = "quay.io/biocontainers/samtools:1.3.1--h1b8c3c0_8"
    //quay.io-biocontainers-samtools-0.1.19--h270b39a_9
  }

  withName:extractZFsfromRawFA{
    cpus = { 1 }
    memory = { 16.GB }
    time = { 2.h * task.attempt}
    container = "docker://kevbrick/prdm9gt_callhaps:1.0"
    //conda = "${params.accessorydir}/conda/callHaps.yml"
  }

  withName:makeHaplotypeTable{
    cpus = { 1 }
    memory = { 8.GB }
    time = { 0.5.h * task.attempt}
    container = "docker://kevbrick/prdm9gt_callhaps:1.0"
  }

  withName:genotypeVsZFlengthDistributions{
    cpus = { 1 }
    memory = { 8.GB }
    time = { 0.5.h * task.attempt}
    container = "docker://kevbrick/prdm9gt_callhaps:1.0"
  }

  withName:drawFigure1{
    cpus = { 1 }
    memory = { 8.GB }
    time = { 0.5.h * task.attempt}
    container = "docker://kevbrick/prdm9gt_callhaps:1.0"
  }

  withName:drawFigure2{
    cpus = { 1 }
    memory = { 8.GB }
    time = { 0.5.h * task.attempt}
    container = "docker://kevbrick/prdm9gt_callhaps:1.0"
  }

  withName:drawFigure3{
    cpus = { 1 }
    memory = { 8.GB }
    time = { 0.5.h * task.attempt}
    container = "docker://kevbrick/prdm9gt_callhaps:1.0"
  }

  withName:drawFigure3MMSupp{
    cpus = { 1 }
    memory = { 8.GB }
    time = { 0.5.h * task.attempt}
    container = "docker://kevbrick/prdm9gt_callhaps:1.0"
  }

  withName:drawGelQuantificationFigure{
    cpus = { 1 }
    memory = { 8.GB }
    time = { 0.5.h * task.attempt}
    container = "docker://kevbrick/gelquant:1.0"
  }

  withName:getMousePrZFAs{
    cpus = { 1 }
    memory = { 8.GB }
    time = { 0.5.h * task.attempt}
    container = "docker://kevbrick/prdm9gt_callhaps:1.0"
  }

  withName:mergePublishedAndFoundAlleles{
    cpus = { 1 }
    memory = { 8.GB }
    time = { 8.h * task.attempt}
    container = "docker://kevbrick/prdm9gt_callhaps:1.0"
  }

  withName:inferRelatednessOfAlleles{
    cpus = { 1 }
    memory = { 8.GB }
    time = { 8.h * task.attempt}
    container = "docker://kevbrick/prdm9gt_callhaps:1.0"
  }

  withName:inferRelatednessOfMouseAlleles{
    cpus = { 1 }
    memory = { 8.GB }
    time = { 6.h * task.attempt}
    container = "docker://kevbrick/prdm9gt_callhaps:1.0"
  }

  withName:processHumanHotspots{
    cpus = { 2 }
    memory = { 8.GB }
    time = { 2.h * task.attempt }
    container = "$baseDir/containers/PrattoEtAl.sif"
  }

  withName:getMotifsChipmunk{
    cpus = { 4 }
    memory = { 8.GB }
    time = { 16.h * task.attempt }
    container = 'kevbrick/chipmunk:1.0'
    conda = ''
  }

  withName:getMotifsMEME{
    cpus = { 4 }
    memory = { 8.GB }
    time = { 2.h * task.attempt }
    container = 'memesuite/memesuite:latest'
  }

}

profiles {
  standard {
    process{
      executor = 'slurm'
      errorStrategy = 'retry'
      maxRetries = 1
      scratch = '/lscratch/$SLURM_JOBID'
      clusterOptions = ' --gres=lscratch:600 '
    }
  }

  local {
    process{
      executor = 'local'
      errorStrategy = 'retry'
      maxRetries = 1
    }
  }
}
