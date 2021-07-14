nextflow.preview.dsl=2

// set some default params
params.help=""

if (params.help) {
  log.info " "
  log.info "=========================================================================="
  log.info "PRDM9 genotyping PIPELINE (Version 4.0)                                "
  log.info "=========================================================================="
  log.info " "
  log.info "USAGE: "
  log.info " "
  log.info "------------------------------------------------------------------------- "
  log.info "nextflow run genotype_prdm9.v4.0.nf "
  log.info " --pb            <string: Pacbio subreads BAMs> "
  log.info " --outdir        <string: default = output> "
  log.info " --pipedir       <string: parent folder for accessoryFiles> "
  log.info " "
  log.info "=========================================================================="
  exit 1
  }

// Params:
params.outdir       = "output"
params.pipedir      = ""
params.inPB         = ""
params.inONT        = ""
params.accessorydir = ""

//log.info
log.info " "
log.info "=========================================================================="
log.info "PRDM9 genotyping PIPELINE (Version 4.0)                                "
log.info "=========================================================================="
log.info " "
log.info "USAGE: "
log.info " "
log.info "------------------------------------------------------------------------- "
log.info "nextflow run genotype_prdm9.nf "
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

process call_ccs {

  publishDir "${params.outdir}/ccs", pattern: '*bam', mode: 'copy', overwrite: false
  publishDir "${params.outdir}/ccs", pattern: '*pbi', mode: 'copy', overwrite: false
  publishDir "${params.outdir}/ccs", pattern: '*txt', mode: 'copy', overwrite: false
  publishDir "${params.outdir}/ccs", pattern: '*tab', mode: 'copy', overwrite: false

  tag { expt }

  input:
  tuple path(rawbam), val(expt), path(f1), path(ccspath), val(dmxpath)

  output:
  tuple(path('*ccs_reads.bam'), path('*ccs_reads.bam.pbi'), val (expt), path('name*tab'), val(dmxpath), emit: bam)
  path('name*tab', emit: nametab)
  path('*ccs_report.txt', emit: report)

  script:
  def outBAM   = expt + ".ccs_reads.bam"
  def outRep   = expt + ".ccs_report.txt"
  """
  ## This allows us to resume from after CCS calling
  if [ -f "${ccspath}/${outBAM}" ]; then
    b=`realpath ${ccspath}/${outBAM}`
    ln -s \$b       ${outBAM}
    ln -s \$b".pbi" ${outBAM}.pbi

    r=`realpath ${ccspath}/${outRep}`
    ln -s \$r       ${outRep}
  else
    ccs ${rawbam} ${outBAM} --report-file ${outRep}
  fi

  grep ${expt} ${f1} >name_conversion_table.tab
  """
  }

process demultiplex {

  //publishDir "${params.outdir}/fa/raw/${expt}", pattern: '*.fa'  , mode: 'copy', overwrite: false
  //publishDir "${params.outdir}/dmx/${expt}",    pattern: '*.bam*', mode: 'copy', overwrite: false
  //publishDir "${params.outdir}/dmx/${expt}",    pattern: '*.xml*', mode: 'copy', overwrite: false
  //publishDir "${params.outdir}/dmx/${expt}",    pattern: '*lima*', mode: 'copy', overwrite: false

  tag { bam }

  input:
  tuple path(bam), path(pbi), val(expt), path(conversiontable), path(dmxpath)
  path(popFile)

  output:
  //path('*raw.fa', emit: fa)
  tuple(path('*bam'),  path('*pbi'),  path('*xml'), val(expt), path(conversiontable), emit: bam)
  path('*lima.counts', emit: dmxCountRep)
  path('*lima.summary', emit: dmxSummary)
  path('all_barcodes_BA.fa', emit: barcodes)

  script:
  def randprefix = new Random().with {(1..30).collect {(('a'..'z')).join()[ nextInt((('a'..'z')).join().length())]}.join()}
  def prefix = bam.name.replaceFirst(/.bam/,"")
  def outDMX = bam.name.replaceFirst(/.bam/,".demux.bam")
  def outRep = bam.name.replaceFirst(/.bam/,".ccs_report.txt")
  """
  echo "${dmxpath}"

  ## This allows us to resume from after demultiplexing
  dmx=`realpath ${dmxpath}`
  fapath=`cd \$dmx; cd ../../fa/raw ; pwd`"/${expt}"

  if [ `find \$dmx -maxdepth 1 -name '*lima.counts' |wc -l` -ge 1 ]; then
    ln -s \$dmx/*bam .
    ln -s \$dmx/*pbi .
    ln -s \$dmx/*xml .
    ln -s \$dmx/*lima* .
    ln -s \$fapath/*fa .
  else
    cp ${params.accessorydir}/barcodes_CCS/* .

    ## Get all combinations of barcodes 1(ONT) & 2(CUSTOM)
    perl ${params.accessorydir}/scripts/combineBCsForCCSDMX.pl >all_barcodes_BA.fa

    lima --ccs --guess 45 --peek 10000 --guess-min-count 5 --different --score-full-pass --split-bam-named ${bam} all_barcodes_BA.fa ${outDMX}

  fi
  """
}

process dmxBAMtoFA {

  publishDir "${params.outdir}/fa/raw/${expt}", pattern: '*.fa'  , mode: 'copy', overwrite: false
  publishDir "${params.outdir}/dmx/${expt}",    pattern: '*.bam*', mode: 'copy', overwrite: false

  tag { bam }

  input:
  tuple path(bam), path(pbi), path(xml), val(expt), path(conversiontable)
  path(popFile)

  output:
  path('*raw.fa', emit: fa)
  tuple(path('*bam'),  path('*pbi'),  path('*xml') , emit: bam)

  shell:
  //  #perl ${params.accessorydir}/scripts/filter_CCS_with_good_demultiplexing.pl --p ${popFile} --n ${conversiontable} --e ${expt}
  '''
  #!/usr/bin/perl
  use strict;

  my $popFile = "!{popFile}";
  my $nameFile = "!{conversiontable}";
  my $exptName = "!{expt}";

  ## GET INDIVIDUAL DETAILS########################################
  my %pop;

  open POP, $popFile;

  while (<POP>){
    chomp;
    my ($pop,$popid,$id,$sex)  = split(/\\t/,$_);

    $pop{$id}->{popid} = $popid;
    $pop{$id}->{pop}   = $pop;
    $pop{$id}->{sex}   = $sex;
  }

  close POP;

  ## GET NAME TO BARCODE TABLE#####################################
  my %names;
  my %edate;

  open NM, $nameFile;

  while (<NM>){
    chomp;
  	next unless ($_ =~ /($exptName)/);

    my ($name,$bc2,$bc1,$expt,$edate)  = split(/,/,$_);

    ## prevent slashes in name
    $name =~ s/[\\/\\\\]/_/g;

    $names{$bc1."_".$bc2} = $name;
    $edate{$expt} = $edate;
  }

  close NM;

  #################################################################
  ## Loop through all demultiplexed BAMs
  opendir(DIR,'.');

  while (my $f = readdir(DIR)){

  	next unless ($f =~ /BC.+BC.+\\.bam$/);

  	my $bam    = $f;
  	my $bamChk = $bam;

  	$bamChk =~ s/BC(\\d+)_([FR])_BC(\\d+).+BC(\\d+)_([FR])_BC(\\d+)//;

  	my ($bc1A,$fr1,$bc1B,$bc2A,$fr2,$bc2B) = ($1,$2,$3,$4,$5,$6);

  	my ($b1,$b2) = ($bc1A."_".$bc1B, $bc2A."_".$bc2B);

  	## DMX good IF: Same names and FR/RF orientation
  	if ($b1 eq $b2 && $fr1 ne $fr2){

  		my $name = $names{$b1};

  		if (!$pop{$name}){
        $pop{$name}->{popid} = 'other';
        $pop{$name}->{pop}   = 'OTH';
        $pop{$name}->{sex}   = 'MALE';
      }

      $name     = join("_",$name,$pop{$name}->{pop},$bc1A,$bc1B,uc($pop{$name}->{sex}) eq "MALE"?"M":"F");

  		my $xml = $bam;
  		$xml =~ s/\\.bam/\\.subreadset.xml/;

  		my $newbam    = $name."_".$exptName.".demux.bam";
  		my $newbampbi = $name."_".$exptName.".demux.bam.pbi";
  		my $newxml    = $name."_".$exptName.".subreadset.xml";
      my $newfa     = $name."_".$exptName.".raw.fa";

  		system("cp $bam $newbam");
  		system("cp $bam.pbi $newbampbi");
     	system("cp $xml ".$name."_".$exptName.".subreadset.xml");

      //system("for f in `samtools view $newbam |cut -f1,10`; do echo $f >>tmp$newfa; done")
      //system("sed \\'1~2s/^/>/\\' $tmpfa |fold -w84 >$newfa")
  		open FA, '>', $newfa;

  		my $seqCnt;
  		open my $IN, '-|', "samtools view $newbam |cut -f1,10";
  		while (<$IN>){
  			chomp;
        my @D = split(/\\t/,$_);
  			print FA ">".$name."_".$exptName."|$D[0]\\n";
  			print FA to84($D[1]);
  		}
  		close FA; close $IN;
  	}else{
  		my $newbam = $bam; $newbam =~ s/demux.bam/demuxBAD.bam/;
  		system("cp $bam $newbam");
  		system("cp $bam.pbi $newbam.pbi");
  		my $xml = $bam;    $xml       =~ s/\\.bam/\\.subreadset.xml/;
  		my $xml = $newbam; my $newxml =~ s/\\.bam/\\.subreadset.xml/;
  		system("cp $xml $newxml");
  	}
  }

  sub to84{
    my $inLine = shift;
    my $lnRet;
    while ($inLine =~ s/^(.{84})//){
      $lnRet .= $1."\\n";
    }
    $lnRet .= $1."\\n";
    return($lnRet);
  }

  '''
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

  for id in `ls *fa |perl -pi -e 's/^(.+).run.*/\$1/' |sort -k1,1 |uniq`; do
    cat \$id*fa >\$id.raw.merged.fa
  done

  ## Split FASTA files to equal sized subfolders
  i=1
  while read l; do
    mkdir -p faSplit_\$i
    ln -s \$l faSplit_\$((i++))
  done< <(find -L `pwd` -name '*merged.fa' | xargs -n 15)
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

//#############################################
// WORKFLOWS
workflow getPacbioRawFastas {
    take: rundataPB

    main:
    prdata    = get_known_human_prdm9_data()

    ccsbam    = call_ccs(rundataPB)
    dmx       = demultiplex(ccsbam.bam, prdata.popData)
    indiv     = dmxBAMtoFA(dmx.bam,prdata.popData)

    emit:
    allfa      = indiv.fa.flatten()
    pubZFs     = prdata.hsZFs
    pubAlleles = prdata.hsAlleles
  }

//#############################################
// MAIN WORKFLOW
workflow {

//  pbIn  = Channel.fromPath(params.inPB).splitCsv(header:false, sep: ",", by: 1)
//                 .map { sample -> tuple(file(sample[0]), sample[1], file(sample[2]), file(sample[3]), file(sample[4])) }

  pbIn  = Channel.fromPath(params.inPB).splitCsv(header:false, sep: ",", by: 1)
                 .map{row ->
                     def input_bam  = row[0]
                     def input_name = row[1]
                     def id_table   = row[2]
                     def ccs_folder = row[3]
                     def dmx_folder = row[4]
                     return [ file(input_bam), input_name, file(id_table), file(ccs_folder), file(dmx_folder) ]}

  raw   = getPacbioRawFastas(pbIn)

  mergeFA     = merge_raw_fa(raw.allfa.collect())
  procFA      = extractZFsfromRawFA(mergeFA.faSplit.flatten(), raw.pubZFs,  raw.pubAlleles)

}
