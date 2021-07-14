nextflow.preview.dsl=2

// set some default params
params.help=""

if (params.help) {
  log.info " "
  log.info "=========================================================================="
  log.info "PRDM9 genotyping PIPELINE (Version 2.0)                                "
  log.info "=========================================================================="
  log.info " "
  log.info "USAGE: "
  log.info " "
  log.info "------------------------------------------------------------------------- "
  log.info "nextflow run genotype_prdm9.nf "
  log.info " --reads         <string: Pacbio subreads BAM> "
  log.info " --outdir        <string: default = output> "
  log.info " --pipedir       <string: parent folder for accessoryFiles> "
  log.info " "
  log.info "=========================================================================="
  exit 1
  }

// Params:
params.outdir         = "output"
params.reads          = ""
params.pipedir        = ""
params.inONT          = ""
params.accessorydir   = ""
params.fast5groupsize = 100
params.mmi            = '/data/RDCO/genomes/hg38/Minimap2/genome.mmi'

def minimap2_index = params.mmi

//log.info
log.info " "
log.info "=========================================================================="
log.info "PRDM9 genotyping PIPELINE (Version 2.0)                                "
log.info "=========================================================================="
log.info " "
log.info "USAGE: "
log.info " "
log.info "------------------------------------------------------------------------- "
log.info "nextflow run genotype_prdm9.nf "
log.info " --reads         ${params.reads} "
log.info " --outdir        ${params.outdir} "
log.info " --pipedir       ${params.pipedir} "
log.info " "
log.info "=========================================================================="

process get_known_human_prdm9_data {

  publishDir "${params.outdir}/publishedData", mode: 'copy', overwrite: true

  output:
  path('humanPRDM9ZFcodes.BergJeffreys.txt',         emit: hsZFs)
  path('humanPRDM9alleles.BergJeffreys.txt',         emit: hsAlleles)
  path('humanPRDM9alleles.BergJeffreys.fa',          emit: hsFA)
  path('humanPRDM9alleles.BergJeffreys.withSeq.txt', emit: hsAlleleSeq)
  path('individuals.1KG.tab',                        emit: popData)

  """
  cp ${params.accessorydir}/scripts/getPRDM9ZFs_JeffreysAndBerg.pl .
  cp ${params.accessorydir}/scripts/replaceUnicode .

  perl getPRDM9ZFs_JeffreysAndBerg.pl

  ## Get 1KG population data
  echo -e "pop\\tpopid\\tid\\tsex" >individuals.1KG.tab

  for pop in "MGP00001","FIN" "MGP00007","TSI" "MGP00008","LWK" "MGP00011","PEL" "MGP00013","YRI" "MGP00017","CHB" "MGP00020","PJL"; do
    IFS=","
    set -- \$pop
    wget https://www.coriell.org/0/Sections/Search/Panel_Detail.aspx?Ref=\$1 -O \$1.tab
    perl -lane 'if (s/^.+((HG|NA)\\d+).+?(female|male).+/"'\$2'\\t'\$1'\\t\$1\\t".lc(\$3)/ie){print \$_}' \$1.tab |sort -k2,2 >>individuals.1KG.tab
  done
  """
  }

process getReadPairs_ONT {

  tag {expt + '_' + f5folder}

  publishDir "${params.outdir}/guppy"    , mode: 'copy', overwrite: true, pattern: 'read_details*txt'

  input:
  tuple path(f5folder), val(expt), path(f1), path(readpairs)
  path(popFile)

  output:
  tuple(path('read_details*txt'), path(f5folder, followLinks: false), val(expt), emit: dets)
  path('*gupp*fastq.gz' ,  optional: true, emit: fq)


  script:
  """
  if [ -f "${readpairs}" ]; then
    cp ${readpairs} read_details.reused.${expt}.txt
  else
    ln -s ${params.accessorydir} .
    f5=`basename ${f5folder}`
    guppy_out="guppy_\$f5"
    guppy_1d2="guppy_1d2_\$f5"

    mkdir \$guppy_out
    mkdir \$guppy_1d2

    ## Get base-calls with guppy
    guppy_basecaller --input_path ${f5folder} -s \$guppy_out --fast5_out -c dna_r9.5_450bps.cfg --device "cuda:all" --disable_pings

    guppy_basecaller_1d2 --compress_fastq -i \$guppy_out/workspace -s \$guppy_1d2 \
                         -c dna_r9.5_450bps_1d2_raw.cfg \
                         --device "cuda:all" --index_file \$guppy_out/sequencing_summary.txt --disable_pings

    ## Parse read pairs detected by guppy
    cut -f3,4 \$guppy_1d2/sequencing_summary.txt |grep -v read_id1 >pairs_by_id.txt

    ## Identify read barcodes (two rounds)
    guppy_barcoder   --compress_fastq -i \$guppy_out -s demux --arrangements_files barcode_arrs_pcr96.cfg --min_score 40 --front_window_size 300 --rear_window_size 300 --trim_barcodes
    echo -e "readname\\tbarcode1\\tbarcode2" >barcodes_by_readname.txt

    for barcode in `ls demux |grep barcode |perl -pi -e 's/barcode//'`; do
        guppy_barcoder --compress_fastq -i demux/barcode\$barcode \
                       -s BC\$barcode \
                       -d accessoryFiles/barcoding_cfgs/ \
                       --arrangements_files custom_12bp.cfg\
                       --min_score 70 --front_window_size 100 --rear_window_size 100 --trim_barcodes

        cut -f1,2 BC\$barcode/barcoding_summary.txt |grep -vP '(unclassified|barcode_arrangement)' |\
             perl -lane '\$bc1="'\$barcode'"; \
                         \$F[1] =~ s/barcode//; print join("\\t",\$F[0],\$bc1,\$F[1])' |sort -k2,2 >>barcodes_by_readname.txt
    done

    grep ${expt} ${f1}  >name_conversion_table.tab

    ## If a symlink already exists, use a different output name
    readDetailsOutput="read_details.${expt}.txt"

    if [ -h "\$readDetailsOutput" ]; then
      unlink \$readDetailsOutput
      readDetailsOutput="read_details.${expt}.new.txt"
    fi

    ## Get a list of "valid" read pairs and associated details (barcodes, indiviudal)
    perl ${params.accessorydir}/scripts/getValidReadPairs.pl --p ${popFile} --n name_conversion_table.tab \
                                                             --bc barcodes_by_readname.txt --pairs pairs_by_id.txt \
                                                             --out \$readDetailsOutput
   fi

   mkdir 1d2Fastq
   find `pwd`/\$guppy_1d2 -name '*fastq.gz' |xargs -I % ln -s % 1d2Fastq

   perl ${params.accessorydir}/scripts/splitFQbyindividual.pl \$readDetailsOutput ${expt}.guppy1d \$guppy_out
   perl ${params.accessorydir}/scripts/splitFQbyindividual.pl \$readDetailsOutput ${expt}.guppy1d2 1d2Fastq

   gzip *fastq ||true
  """
  }

process getBarcodes_ONT {

  tag {expt + '_' + f5folder}

  publishDir "${params.outdir}/guppy"    , mode: 'copy', overwrite: true
  input:
  tuple path(f5folder), val(expt), path(f1), path(readpairs)
  path(popFile)

  output:
  //tuple(path('read_details*txt'), path(f5folder, followLinks: false), val(expt), emit: dets)
  tuple(val(expt), path('reads_by_individual.txt'), path(f5folder, followLinks: false), emit: dets)

  script:
  """
  if [ -f "${readpairs}" ]; then
    cp ${readpairs} read_details.reused.${expt}.txt
  else
    ln -s ${params.accessorydir} .
    f5=`basename ${f5folder}`
    guppy_out="guppy_\$f5"

    mkdir \$guppy_out

    ## Get base-calls with guppy
    guppy_basecaller --input_path ${f5folder} -s \$guppy_out --fast5_out -c dna_r9.5_450bps.cfg --device "cuda:all" --disable_pings

    ## Identify read barcodes (two rounds)
    guppy_barcoder   --compress_fastq -i \$guppy_out -s demux --arrangements_files barcode_arrs_pcr96.cfg --min_score 40 --front_window_size 300 --rear_window_size 300 --trim_barcodes
    echo -e "readname\\tbarcode1\\tbarcode2" >barcodes_by_readname.txt

    for barcode in `ls demux |grep barcode |perl -pi -e 's/barcode//'`; do
        guppy_barcoder --compress_fastq -i demux/barcode\$barcode \
                       -s BC\$barcode \
                       -d accessoryFiles/barcoding_cfgs/ \
                       --arrangements_files custom_12bp.cfg\
                       --min_score 70 --front_window_size 100 --rear_window_size 100 --trim_barcodes

        cut -f1,2 BC\$barcode/barcoding_summary.txt |grep -vP '(unclassified|barcode_arrangement)' |\
             perl -lane '\$bc1="'\$barcode'"; \
                         \$F[1] =~ s/barcode//; print join("\\t",\$F[0],\$bc1,\$F[1])' |sort -k2,2 >>barcodes_by_readname.txt
    done

    cat barcodes_by_readname.txt |perl -lane 'print join("\\t",\$F[1].":".\$F[2],\$F[0])'                                   |sort -k1,1 >barcodes_reads.txt
    grep ${expt} ${f1} |perl -F, -lane '\$F[2] = "0\$F[2]" if (length(\$F[2]) < 2); print join("\\t",\$F[2].":".\$F[1],\$F[0])' |sort -k1,1 >name_conversion_table.tab
    join -1 1 -2 1 barcodes_reads.txt name_conversion_table.tab |awk -F" " '{print \$3,\$2}' OFS="\\t" >reads_by_individual.txt
  fi
  """
  }

process splitIndividualsForBonito {

    tag {readPairDetails}

    input:
    tuple path(readPairDetails), path(f5folder), val(expt)

    output:
    path('read_dets_*txt', emit: readList)
    tuple(val(expt), path(f5folder, followLinks: false), emit: all)

    script:
    """
    grep -v UNK  ${readPairDetails} |grep -v individual |sort -k7,7 |cut -f7 |uniq |split -l ${params.fast5groupsize} - read_dets_

    for r in read_dets_*; do
      grep -Ff \$r ${readPairDetails} |sort -k7,7 >\$r"_${expt}.txt"
      ln -s `realpath ${f5folder}` \$r"_${expt}_fast5"
    done
    """
  }

process bonitoPairDecoding_ONT {

  tag {expt}

  publishDir "${params.outdir}/fasta/bonito/raw", mode: 'copy', overwrite: true, pattern: '*fa'

  input:
  tuple  val(expt), path(readPairDetails), path(f5folder)

  output:
  path('*.raw.fa', emit: fa)

  script:
  """
  ln -s ${params.accessorydir} .

  ## Make read-pairs files
  perl -lane 'print \$F[0]." ".\$F[1] if (\$_ !~ /(name|UNK|unk)/);'   ${readPairDetails} >read_pairs.txt
  perl -lane 'print \$F[0]."\\n".\$F[1] if (\$_ !~ /(name|UNK|unk)/);' ${readPairDetails} >allreads.txt

  ## Make new fast5s with ONLY the reads for 1D2
  ## Takes a few minutes but saves TONS of time on basecalling
  mkdir subset
  fast5_subset -i ${f5folder} -t 6 -s subset -f ${expt}_validPairs_ -n 1000 -l allreads.txt

  ## Do base-calling
  bonito pair --half read_pairs.txt subset >all.fa

  ## extract per-individual sequences from fasta using a list of IDs
  for nm in `cut -f7 ${readPairDetails} |grep -vP '(unk|UNK|name)' |sort |uniq`; do

    grep \$nm ${readPairDetails} |cut -f1,2 |perl -lane 'print \$F[0].";".\$F[1].";"' >\$nm"_ids.txt"

    awk -F'>' 'NR==FNR{ids[\$0]; next} NF>1{f=(\$2 in ids)} f' \$nm"_ids.txt" all.fa >\$nm"_${expt}.raw.fa"

    ## This may make some empty files. If so, delete them ...
    if [ -s \$nm"_${expt}.raw.fa" ]; then
      echo \$nm"_${expt}.raw.fa ... OK"
    else
      rm \$nm"_${expt}.raw.fa"
    fi

  done

  """
  }

process bonitoBasecallSR_ONT {

  tag {expt}

  publishDir "${params.outdir}/fasta/bonito/raw", mode: 'copy', overwrite: true, pattern: '*fa'

  input:
  tuple  val(expt), path(readPairDetails), path(f5folder)

  output:
  path('*.raw.fa', emit: fa)

  script:
  """
  ln -s ${params.accessorydir} .

  ## Make read-pairs files
  cut -f2 ${readPairDetails} >allreads.txt

  ## Do base-calling
  bonito basecaller dna_r9.4.1 --read-ids allreads.txt --device cuda:0 `pwd`/${expt} |paste -d" " - - |cut -c2- >basecalls.txt

  ## Add individual name & split to individual fasta files
  sort -k2,2 ${readPairDetails} >rbi_sorted.txt
  sort -k1,1 basecalls.txt >basecalls.sorted.txt
  join -1 2 -2 1 rbi_sorted.txt basecalls.sorted.txt |perl -pi -e 's/\\//-/' |awk '{print>\$2".basecalls.txt"}' 
  
  for sam in *basecalls.txt; do 
	fa=\${sam/basecalls.txt/raw.fa}
        cat \$sam | perl -lane 'print ">\$F[0]|\$F[1]\\n\$F[2]"' |fold -w 84  >\$fa
  done

  """
  }

process merge_raw_fa {
  publishDir "${params.outdir}/fasta/rawmerged", pattern: '*.fa', mode: 'copy', overwrite: true

  input:
  path(fa)

  output:
  path('*.fa'     , emit: fa)
  path('faSplit*' , emit: faSplit)

  """
  ln -s ${params.accessorydir} .

  for id in `ls *fa |perl -pi -e 's/^(\\S+?_\\S\\S\\S)_.+\$/\$1/' |sort -k1,1 |uniq`; do
    cat \$id*fa >\$id.raw.merged.fa
  done

  ## Split FASTA files to equal sized subfolders
  i=1
  while read l; do
    mkdir -p faSplit_\$i
    ln -s \$l faSplit_\$((i++))
  done< <(find -L `pwd` -name '*merged.fa' | xargs -n 30)
  """
  }

process guppyFQtoFA {
  publishDir "${params.outdir}/fasta/rawguppy/1D" , pattern: '*guppy1d.merge.raw.fa'     , mode: 'copy', overwrite: true
  publishDir "${params.outdir}/fasta/rawguppy/1D2", pattern: '*guppy1d2.merge.raw.fa'    , mode: 'copy', overwrite: true
  publishDir "${params.outdir}/fastq/rawguppy/1D" , pattern: '*guppy1d.merge.fastq.gz'   , mode: 'copy', overwrite: true
  publishDir "${params.outdir}/fastq/rawguppy/1D2", pattern: '*guppy1d2.merge.fastq.gz'  , mode: 'copy', overwrite: true

  input:
  path(fq)

  output:
  path('*.fa'        , emit: fa)
  path('*.fastq.gz'  , emit: fq)

  """
  for individual in `ls *1d*gz |perl -pi -e 's/^(.+?)\\..+\$/\$1/' |sort |uniq`; do
    if [ `ls \$individual*guppy1d.fastq.gz |wc -l` -gt 0 ]; then
      zcat \$individual*guppy1d.fastq.gz   |gzip -c >\$individual".guppy1d.merge.fastq.gz"
      zcat \$individual".guppy1d.merge.fastq.gz" |paste - - - - | cut -f 1,2 | sed 's/^@/>/' | tr "\\t" "\\n" > \$individual".guppy1d.merge.raw.fa"
    fi

    if [ `ls \$individual*guppy1d2.fastq.gz |wc -l` -gt 0 ]; then
      zcat \$individual*guppy1d2.fastq.gz   |gzip -c >\$individual".guppy1d2.merge.fastq.gz"
      zcat \$individual".guppy1d2.merge.fastq.gz" |paste - - - - | cut -f 1,2 | sed 's/^@/>/' | tr "\\t" "\\n" > \$individual".guppy1d2.merge.raw.fa"
    fi
  done

  """
  }

process extractZFsfromRawFA{
  publishDir "${params.outdir}/fasta/ZF",   pattern: '*.fa',     mode: 'copy', overwrite: true
  publishDir "${params.outdir}/haplotypes", pattern: '*haplot*', mode: 'copy', overwrite: true
  publishDir "${params.outdir}/debug",      pattern: '*scores*', mode: 'copy', overwrite: true
  publishDir "${params.outdir}/debug",      pattern: '*errorp*', mode: 'copy', overwrite: true

  input:
  path(faPath)
  path(pubZFs)
  path(pubAlleles)

  output:
  path('*ZFs.fa',               emit: fa)
  path('*haplotypes.txt',       emit: hap)
  path('*consensus_scores.txt', emit: sc)
  path('*errorprofile.txt',     emit: err)

  """
  ## Assure we have BioPerl on perl5lib path
  export PERL5LIB=\$PERL5LIB:\$CONDA_PREFIX"/lib/perl5/site_perl/5.22.0/"

  for fa in ${faPath}/*.fa; do
    fnm=`basename \$fa`
    bcid=\${fnm/.raw.merged.fa/}
    id=\$bcid

    perl ${params.accessorydir}/scripts/extractPRDM9ZFsFromFA.pl \
      --fa \$fa --id \$id --rspath ${params.accessorydir}/scripts \
      --type \$id --n 30 \
      --ZFs ${pubZFs} --aln blast

  	nfa=`ls \$id*ZFs.fa | wc -l`
  	onehap=''

  	if [ \$nfa == 2 ]; then
  		onehap=' --onehap '
  	fi

  	if [[ \$nfa -ge 1 ]]; then
  		for zfa in \$id*ZFs.fa; do
  			perl ${params.accessorydir}/scripts/getZFAhaplotypesFromFA.pl --fa \$zfa \$onehap --id \$id --pa ${pubAlleles} --pz ${pubZFs}
  		done
  	else
  		echo -e "\$bcid\\tnoZFA\\tNA\\tNA\\tNA\\tNA" >\$id.haplotypes.txt
  	fi

    mkdir -p zfs_\$id
    cp \$bcid*.* zfs_\$id

  done
  """
  }

// OK ... let's start
workflow getONTRawFastas1D2{
  take: rundata

  main:
  prdata      = get_known_human_prdm9_data()

  readpairs   = getReadPairs_ONT(rundata,prdata.popData)
  initSplit   = splitIndividualsForBonito(readpairs.dets)

  f5Splits = initSplit.readList.flatten().map { sample -> tuple(sample.name.replaceFirst('^read_dets_[a-z]{2}_(.+)(.txt)$','\$1'),
                                                                file(sample))}
                         .combine(initSplit.all, by: 0)

  guppyRaw    = guppyFQtoFA(readpairs.fq.collect())
  bonito      = bonitoPairDecoding_ONT(f5Splits)

  emit:
  allfa      = bonito.fa.flatten()
  pubZFs     = prdata.hsZFs
  pubAlleles = prdata.hsAlleles
  }

// OK ... let's start
workflow getONTRawFastas{
  take: rundata

  main:
  prdata      = get_known_human_prdm9_data()

  readpairs   = getBarcodes_ONT(rundata,prdata.popData)

  bonito      = bonitoBasecallSR_ONT(readpairs)

  emit:
  allfa      = bonito.fa.flatten()
  pubZFs     = prdata.hsZFs
  pubAlleles = prdata.hsAlleles
  }

//#############################################
// MAIN WORKFLOW
workflow {

  ontIn  = Channel.fromPath(params.inONT).splitCsv(header:false, sep: ",", by: 1)
                 .map { sample -> tuple(file(sample[0]), sample[1], file(sample[2]), file(sample[3])) }

  ontIn.view()

  raw         = getONTRawFastas(ontIn)
}
