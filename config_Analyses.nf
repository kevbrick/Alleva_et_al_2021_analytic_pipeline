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
singularity.runOptions=" -B \$NXF_GENOMES -B ${params.accessorydir} -B ${params.pipedir} -B ${launchDir}"

process{

  withName:get_known_human_prdm9_data{
    cpus = { 2 }
    memory = { 8.GB }
    time = { 0.5.h * task.attempt }
    //container = "docker://kevbrick/parsetools:1.0"
  }

  withName:merge_raw_fa{
    cpus = { 1 }
    memory = { 8.GB }
    time = { 2.h * task.attempt}
  }

  withName:extractZFsfromRawFA{
    cpus = { 1 }
    memory = { 16.GB }
    time = { 9.h * task.attempt}
    container = "docker://kevbrick/prdm9gt_callhaps:1.0"
  }

  withName:makeHaplotypeTable{
    cpus = { 1 }
    memory = { 8.GB }
    time = { 0.5.h * task.attempt}
    container = "docker://kevbrick/prdm9gt_callhaps:1.0"
  }

  withName:checkTrios{
    cpus = { 2 }
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
  
  withName:parseNewGTsONLY{
    cpus = { 1 }
    memory = { 2.GB }
    time = { 0.25.h * task.attempt}
    container = "docker://kevbrick/alleva_shortreadvalidation:1.0"
  }

  withName:checkNewAllelesWithShortReads{
    cpus = { 4 }
    memory = { 16.GB }
    time = { 4.h * task.attempt}
    container = "docker://kevbrick/alleva_shortreadvalidation:1.0"
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
    container = "docker://kevbrick/alleva_r:2.0"
  }

  withName:drawRelatednessFigure{
    cpus = { 1 }
    memory = { 8.GB }
    time = { 1.h * task.attempt}
    container = "docker://kevbrick/prdm9gt_callhaps:1.0"
  }

  withName:drawShortReadValidationFig{
    cpus = { 1 }
    memory = { 8.GB }
    time = { 1.h * task.attempt}
    container = "docker://kevbrick/prdm9gt_callhaps:1.0"
  }
  
  withName:drawHotspotsFigure{
    cpus = { 1 }
    memory = { 8.GB }
    time = { 0.5.h * task.attempt}
    container = "docker://kevbrick/alleva_r:2.0"
  }

  withName:drawGelQuantificationFigure{
    cpus = { 1 }
    memory = { 8.GB }
    time = { 0.5.h * task.attempt}
    container = "docker://kevbrick/gelquant:1.0"
  }

  withName:analyzeDiscordancesPB_v_ONT{
    cpus = { 1 }
    memory = { 8.GB }
    time = { 1.h * task.attempt}
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
    time = { 28.h * task.attempt}
    container = "docker://kevbrick/alleva_inferrelatedness:1.0"
  }

  withName:getLinkedAlleles{
    cpus = { 4 }
    memory = { 16.GB }
    time = { 2.h * task.attempt}
    container = "docker://kevbrick/alleva_plink:2.0"
  }

  withName:drawAssociationsFigure{
    cpus = { 1 }
    memory = { 8.GB }
    time = { 1.h * task.attempt}
    container = "docker://kevbrick/prdm9gt_callhaps:1.0"
  }

  withName:dnaToPeptide{
    cpus = { 1 }
    memory = { 8.GB }
    time = { 1.h * task.attempt}
    container = "docker://kevbrick/prdm9gt_callhaps:1.0"
  }

  withName:makeACtypesPlot{
    cpus = { 1 }
    memory = { 8.GB }
    time = { 2.h * task.attempt}
    container = "docker://kevbrick/alleva_r:2.0"
  }

  withName:processHumanHotspots{
    cpus = { 2 }
    memory = { 8.GB }
    time = { 2.h * task.attempt }
    container = "docker://kevbrick/alleva_r:2.0"
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
