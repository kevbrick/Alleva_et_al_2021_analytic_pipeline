params.pipedir="."
params.accessorydir = "${params.pipedir}/accessoryFiles/"

timeline {
  enabled = true
  file = "${params.outdir}/nextflow/prdm9_ont_1kg_timeline.html"
}

report {
  enabled = true
  file = "${params.outdir}/nextflow/prdm9_ont_1kg_report.html"
}

trace {
  enabled = true
  file = "${params.outdir}/nextflow/prdm9_ont_1kg_trace.txt"
}

manifest {
  description = 'Human 1KG PRDM9 genotyping pipeline (ONT). Author: Kevin Brick.'
}

singularity.enabled = true
singularity.autoMounts = true
singularity.envWhitelist='https_proxy,http_proxy,ftp_proxy,DISPLAY,NXF_GENOMES'
singularity.runOptions=" -B \$NXF_GENOMES -B ${params.accessorydir} -B ${params.pipedir} --nv"

process{

  withName:getReadPairs_ONT{
    time           = { 24.h * task.attempt }
    container      = "docker://kevbrick/ont_guppybonito:3.1.5_0.23.0"
    //module = 'guppy/3.1.5'
    clusterOptions = ' --partition=gpu --cpus-per-task=6 --mem=16g --gres=lscratch:500,gpu:p100:1 '
  }

  withName:splitIndividualsForBonito{
    time           = { 3.h * task.attempt }
    container      = "docker://kevbrick/ont_guppybonito:3.1.5_0.23.0"
    clusterOptions = ' --partition=gpu --cpus-per-task=6 --mem=16g --gres=lscratch:500,gpu:p100:1 '
  }

  withName:splitFAST5sByIndividual{
    time           = { 8.h * task.attempt }
    container      = "docker://kevbrick/ont_guppybonito:3.1.5_0.23.0"
    clusterOptions = ' --partition=gpu --cpus-per-task=6 --mem=16g --gres=lscratch:500,gpu:p100:1 '
  }

  withName:bonitoBasecall_ONT{
    time           = { 16.h * task.attempt }
    container      = "docker://kevbrick/ont_guppybonito:3.1.5_0.23.0"
    clusterOptions = ' --partition=gpu --cpus-per-task=6 --mem=16g --gres=lscratch:500,gpu:p100:1 '
  }

  withName:split_ONT1d_by_barcode{
    cpus = { 2 }
    memory = { 8.GB }
    time = { 0.5.h * task.attempt }
  }

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

}

profiles {
  standard {
    process{
      executor = 'slurm'
      errorStrategy = 'retry'
      maxRetries = 1
      scratch = '/lscratch/$SLURM_JOBID'
      clusterOptions = ' --gres=lscratch:600 '
      pollInterval = '1 min'
      queueStatInterval = '2 min'
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
