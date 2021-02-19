params.pipedir="."
//params.outdir="${launchDir}/${params.outdir}"
params.accessorydir = "${params.pipedir}/accessoryFiles/"

timeline {
  enabled = false
  file = "${params.outdir}/nextflow/prdm9_genotyping_1kg_timeline.html"
}

report {
  enabled = false
  file = "${params.outdir}/nextflow/prdm9_genotyping_1kg_report.html"
}

trace {
  enabled = false
  file = "${params.outdir}/nextflow/prdm9_genotyping_1kg_trace.txt"
}

manifest {
  description = 'Human 1KG PRDM9 genotyping pipeline. Author: Kevin Brick.'
}

singularity.enabled = true
singularity.autoMounts = true
singularity.envWhitelist='https_proxy,http_proxy,ftp_proxy,DISPLAY,NXF_GENOMES'
singularity.runOptions=" -B \$NXF_GENOMES -B ${params.accessorydir} -B ${params.pipedir} -B ${launchDir}/${params.outdir}"

process{

  withName:get_known_human_prdm9_data{
    cpus = { 2 }
    memory = { 8.GB }
    time = { 0.5.h * task.attempt }
  }

  withName:merge_raw_fa{
    cpus = { 1 }
    memory = { 8.GB }
    time = { 2.h * task.attempt}
  }

  withName:extractZFsfromRawFA{
    cpus = { 1 }
    memory = { 16.GB }
    time = { 4.h * task.attempt}
    container = "docker://kevbrick/prdm9gt_callhaps:1.0"
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
    memory = { 16.GB }
    time = { 4.h * task.attempt}
    container = "docker://kevbrick/prdm9gt_callhaps:1.0"
  }

  withName:drawFigure4{
    cpus = { 1 }
    memory = { 8.GB }
    time = { 0.5.h * task.attempt}
    container = "docker://kevbrick/alleva_r:1.0"
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
    container = "docker://kevbrick/alleva_r:1.0"
  }

  withName:getMotifsChipmunk{
    cpus = { 4 }
    memory = { 8.GB }
    time = { 16.h * task.attempt }
    container = 'docker://kevbrick//chipmunk:1.0'
  }

  withName:getMotifsMEME{
    cpus = { 4 }
    memory = { 8.GB }
    time = { 2.h * task.attempt }
    container = 'docker://memesuite/memesuite:latest'
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
