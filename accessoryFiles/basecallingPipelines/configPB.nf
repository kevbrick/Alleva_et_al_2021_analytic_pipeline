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

process{

  withName:getNamedExptChannelPB{
    cpus = { 2 }
    memory = { 8.GB }
    time = { 0.5.h * task.attempt }
    process.conda = "${params.accessorydir}/conda/environment.yml"
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
    process.conda = "${params.accessorydir}/conda/environment.yml"
  }

  withName:call_ccs{
    cpus = { 16 }
    memory = { 32.GB }
    time = { 16.h * task.attempt }
    process.conda = "${params.accessorydir}/conda/environment.yml"
  }

  withName:demultiplex{
    cpus = { 16 }
    memory = { 32.GB }
    time = { 16.h * task.attempt}
    process.conda = "${params.accessorydir}/conda/environment.yml"
  }

  withName:gatherCCSDemuxFAs{
    cpus = { 2 }
    memory = { 8.GB }
    time = { 3.h * task.attempt}
    process.conda = "${params.accessorydir}/conda/environment.yml"
  }

  withName:merge_raw_fa{
    cpus = { 1 }
    memory = { 8.GB }
    time = { 2.h * task.attempt}
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
