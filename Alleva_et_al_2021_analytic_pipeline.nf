nextflow.preview.dsl=2

// set some default params
params.help=""

if (params.help) {
  log.info " "
  log.info "=========================================================================="
  log.info "PRDM9 genotyping PIPELINE (Version 1.0)                                "
  log.info "=========================================================================="
  log.info " "
  log.info "USAGE: "
  log.info " "
  log.info "------------------------------------------------------------------------- "
  log.info "nextflow run Alleva_et_al_2021_analytic_pipeline.nf "
  log.info " --reads         <string: Pacbio subreads BAM> "
  log.info " --outdir        <string: default = output> "
  log.info " --pipedir       <string: parent folder for accessoryFiles> "
  log.info " "
  log.info "=========================================================================="
  exit 1
  }

// Params:
params.outdir         = ""
params.pbfa           = ""
params.ontfa          = ""
params.accessorydir   = ""

params.genomefa       = ""
params.genomeidx      = ""
params.npeaks         = 1000
params.bigwinSize     = 2000
params.centerSz       = 1000
params.RM             = true

genomeFA  = file("${params.genomefa}")
genomeIDX = file("${params.genomefa}.fai")

//log.info
log.info " "
log.info "=========================================================================="
log.info "PRDM9 genotyping PIPELINE for Alleva et al. 2021 (Version 2.0)            "
log.info "=========================================================================="
log.info " "
log.info "USAGE: "
log.info " "
log.info "------------------------------------------------------------------------- "
log.info "nextflow run Alleva_et_al_2021_analytic_pipeline.nf "
log.info " --pbfa             ${params.pbfa} "
log.info " --ontfa            ${params.ontfa} "
log.info " --outdir           ${params.outdir} "
log.info " --accessorydir     ${params.accessorydir} "
log.info " "
log.info "=========================================================================="

process get_known_human_prdm9_data {

  publishDir "output/publishedData", mode: 'copy', overwrite: true

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

process get_jeffreys_blood_and_sperm_variants_only {

  publishDir "output/publishedData", mode: 'copy', overwrite: true

  output:
  path('humanPRDM9alleles.BloodandSpermVariants.txt',  emit: alleles)

  """
  cp ${params.accessorydir}/scripts/getPRDM9ZFs_JeffreysAndBerg.pl .
  cp ${params.accessorydir}/scripts/replaceUnicode .

  perl getPRDM9ZFs_JeffreysAndBerg.pl
  """
  }

process merge_raw_fa {
  publishDir "${params.outdir}/fasta/rawmerged", pattern: '*.fa', mode: 'copy', overwrite: true

  input:
  path(fa)
  val(type)

  output:
  path('*raw.fa', emit: fa)
  //path('faSplit*' ,       emit: faSplit)

  """
  ln -s ${params.accessorydir} .

  for id in `ls *fa |perl -pi -e 's/^(.+?_\\S\\S\\S).+\$/\$1/' |sort -k1,1 |uniq`; do

    cat \$id*fa >\$id.${type}.raw.fa

  done
  """

  }

process extractZFsfromRawFA{
  publishDir "${params.outdir}/fasta/${type}/ZF",   pattern: '*.fa',         mode: 'copy', overwrite: true
  publishDir "${params.outdir}/${type}/haplotypes", pattern: '*haplot*',     mode: 'copy', overwrite: true
  publishDir "${params.outdir}/${type}/debug",      pattern: '*.*',          mode: 'copy', overwrite: true
  publishDir "${params.outdir}/${type}/dets",       pattern: '*zf_arr*det*', mode: 'copy', overwrite: true

  input:
  val(type)
  path(faPath)
  path(pubZFs)
  path(pubAlleles)

  output:
  path('*ZFs.fa',               optional: true, emit: fa)
  path('*haplotypes.txt',       optional: true, emit: hap)
  path('*consensus_scores.txt', optional: true, emit: sc)
  path('*errorprofile.txt',     optional: true, emit: err)
  path('*zf_arr*details.txt',   optional: true, emit: dets)
  path('zfs_*',                 optional: true, emit: debug)

  """
  ## Assure we have BioPerl on perl5lib path
  ## export PERL5LIB=\$PERL5LIB:\$CONDA_PREFIX"/lib/perl5/site_perl/5.22.0/"
  export PERL5LIB="/opt/conda/envs/parsePRDM9RE/lib/site_perl/5.26.2"

  for fa in *.fa; do 
    id=`basename \$fa |perl -pi -e 's/^(\\S+?_\\S\\S\\S).+\$/\$1/'`

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
  		echo -e "\$id\\tnoZFA\\tNA\\tNA\\tNA\\tNA" >\$id.haplotypes.txt
  	fi

    mkdir -p zfs_\$id

    find -maxdepth 1 -name '\$id*.*' -exec cp {} ./zfs_\$id/ \\;

  done
  """
  }

process makeHaplotypeTable{
  publishDir "${params.outdir}/haplotypesTable/${type}", mode: 'copy', overwrite: true

  input:
  val(type)
  path(haps)
  path(fa)
  path(zffa)

  output:
  path('prdm9_haplotypes.*.tab', emit: hap)
  path('new*',                   emit: newZF)

  """
  perl ${params.accessorydir}/scripts/gatherHaps.pl >allhap.txt

  echo -e `grep allele_1_code allhap.txt`"\\tnseqs\\tnzfseqs" |perl -pi -e 's/\\s+(\\S)/\\t\$1/g' >prdm9_haplotypes.${type}.tmp

  for f in *.fa; do

    if [[ ! \$f =~ "ZFs.fa" ]]; then
      nm=`basename \$f |perl -pi -e 's/^(\\S+?)_.+\$/\$1/' 2>/dev/null`

      if [[ `find -name "\$nm*ZFs.fa"` ]]; then
        sz=`grep -P "^\\>" \$f |wc -l`
        zfsz=`cat \$nm*ZFs.fa |grep \\> |wc -l`

        hap=`grep \$nm allhap.txt |head -n1`

        echo -e "\$hap\\t\$sz\\t\$zfsz" >>prdm9_haplotypes.${type}.tmp
        echo "Processed \$nm ..."
      else
        echo "Skipping \$nm: No ZFs found ..."
      fi
    fi
  done

  grep YRI ${params.accessorydir}/pops/20140610_all_samples_with_children_no_quotes.ped |
         perl -lane 'BEGIN{print join("\\t","child","mom","dad")};
                     print join("\\t",@F[1..3]) if (\$F[2] && \$F[3])' >trio.txt

  grep YRI ${params.accessorydir}/pops/20140610_all_samples_with_children_no_quotes.ped |
        perl -lane 'print \$F[1] if (\$F[2] && \$F[3])' >kids.txt

  grep -v -f kids.txt prdm9_haplotypes.${type}.tmp |perl -lane 'push @F, (\$_ =~ /allele_1_code/)?"is_child":"FALSE"; print join("\\t",@F)'  >withkids.tmp
  grep    -f kids.txt prdm9_haplotypes.${type}.tmp |perl -lane 'push @F, "TRUE"; print join("\\t",@F)'                                      >>withkids.tmp

  grep allele_1_code withkids.tmp                         >prdm9_haplotypes.${type}.tab
  sort -k1,1         withkids.tmp |grep -v allele_1_code >>prdm9_haplotypes.${type}.tab

  touch newAlleles.txt
  touch newZFs.txt
  """
  }

process checkTrios{
  publishDir "${params.outdir}/haplotypesTable", mode: 'copy', overwrite: true

  input:
  path(haps)

  output:
  path('trio_haplotypes.tab', emit: hap)

  """
  ln -s ${params.accessorydir} accessoryFiles

  cp accessoryFiles/otherdata/mouse_PrZFA_alleles_in_Pubs.txt .
  
  grep YRI accessoryFiles/pops/20140610_all_samples_with_children_no_quotes.ped |
         perl -lane 'BEGIN{print join("\\t","child","mom","dad")};
                     print join("\\t",@F[1..3]) if (\$F[2] && \$F[3])' >trio.txt
  
  perl accessoryFiles/scripts/checkTrios.pl ${haps} trio.txt >trio_haplotypes.tab
  
  """
  }
  
process genotypeVsZFlengthDistributions{
  publishDir "${params.outdir}/ZFlengths", mode: 'copy', overwrite: true

  input:
  val(type)
  path(haps_table)
  path(dets)

  output:
  path('allZFlengths.*.tab', emit: tab)

  """
  echo -e "id\tpop\tzf_len\tfreq" >allZFlengths.${type}.tab

  for f in *_zf_ar*details.txt; do
    b=`basename \$f`
    nm=`echo "\$b" | cut -d'.' -f1`
    perl -lane 'if (\$_ =~ /ZF\\sCOUNT:\\s+(\\d+)/){print \$1}' \$f |sort -k1n,1n |uniq -c |perl -pi -e 's/\\s+(\\d+)\\s+(\\d+)/'\$nm'\\t\$2\\t\$1/' |perl -pi -e 's/^(\\S+)_(\\S+)/\$1\\t\$2/' 2>/dev/null >>allZFlengths.${type}.tab
  done
  """
  }

process analyzeDiscordancesPB_v_ONT{
  publishDir "${params.outdir}/figures",    mode: 'copy', overwrite: true, pattern: '*png'
  publishDir "${params.outdir}/figures",    mode: 'copy', overwrite: true, pattern: '*pdf'

  input:
  path(ont_dets)
  path(pacbio_dets)

  output:
  path('Allev*.p??', emit: mmFig)

  """
  mkdir nanopore
  mkdir pacbio

  ln -s ${launchDir}/${params.outdir}/ont/dets/*details.txt ./nanopore
  ln -s ${launchDir}/${params.outdir}/pacbio/dets/*details.txt ./pacbio

  ## GET ZFA sizes
  echo -e "id\\tpop\\tsize" >nanopore_ZFAsizes.txt
  echo -e "id\\tpop\\tsize" >pacbio_ZFAsizes.txt

  for x in nanopore/*.prdm9_zf_array*details.txt; do
    nm=`basename \$x`
    name=`echo \${nm/.prdm9_zf_array*/} | perl -pi -e 's/_/","/' 2>/dev/null`
    cat \$x |perl -lane 'if (\$_ =~ /ZF\\sCOUNT:\\s+(\\d+)/){
                           print join("\\t","'\$name'",\$1)}' >>nanopore_ZFAsizes.txt
  done

  for x in pacbio/*.prdm9_zf_array*details.txt; do
    nm=`basename \$x`
    name=`echo \${nm/.prdm9_zf_array*/} | perl -pi -e 's/_/","/' 2>/dev/null`
    cat \$x |perl -lane 'if (\$_ =~ /ZF\\sCOUNT:\\s+(\\d+)/){
                           print join("\\t","'\$name'",\$1)}' >>pacbio_ZFAsizes.txt
  done

  R --no-save <${params.accessorydir}/scripts/drawDiscordantCallsFig.R
  rm -f Rplots.pdf
  """
  }
 
process parseNewGTsONLY {
  input:
  path(haplotypes)

  output:
  path("individuals_with_new_alleles.tab", emit: tab)

  script:
  """
  grep -P '\\sM\\d+\\s' ${haplotypes} >individuals_with_new_alleles.tab
  """
  }

process checkNewAllelesWithShortReads {

  input:
  path(haplotypes)
  path(dets)

  output:
  path("shortreadvalidation*tab", emit: tab)

  script:
  def randID = new Random().with {(1..30).collect {(('a'..'z')).join()[ nextInt((('a'..'z')).join().length())]}.join()}
  """
  ln -s ${params.accessorydir} accessoryFiles

  cp accessoryFiles/scripts/validateNewAllelesWithShortReadData.pl .
  #cp accessoryFiles/1KG/ftp_1KG.current.tree.June232021.txt ./ftp.current.tree
  cp accessoryFiles/1KG/gs_1KG.tree ./ftp.current.tree
  
  perl validateNewAllelesWithShortReadData.pl --haps ${haplotypes} \
                   --noprefix \
                   --tree ftp.current.tree \
                   --dets ${dets} >shortreadvalidation.${randID}.tab 
  """
  }

process drawShortReadValidationFig {

  publishDir "${params.outdir}/figures", mode: 'copy', overwrite: true
  
  input:
  path(dets)

  output:
  path("Alleva*p??", emit: img)

  script:
  """
  ln -s ${params.accessorydir} accessoryFiles

  cat shortreadvalidation*tab >allShortReadValidationData.tab
  
  R --no-save <${params.accessorydir}/scripts/drawAllevaShortReadValidationPlots.R
  
  rm -f Rplots.pdf
  """
  }

process drawFigure1{
  publishDir "${params.outdir}/figures", mode: 'copy', overwrite: true

  input:
  path(ont_dets)
  path(pacbio_dets)

  output:
  path('*png', emit: png)
  path('*pdf', emit: pdf)

  """
  ln -s ${params.accessorydir} accessoryFiles 
  
  cp prdm9_haplotypes.ont.tab prdm9_haplotypes.bonito.tab
  
  R --no-save <${params.accessorydir}/scripts/drawAllevaFigure1.R
  rm -f Rplots.pdf
  """
  }

process drawFigure2{
  publishDir "${params.outdir}/figures", mode: 'copy', overwrite: true

  input:
  path(dets)
  path(allele_AorCtypes)

  output:
  path('Alleva*png', emit: png)
  path('Alleva*pdf', emit: pdf)

  """
  ln -s ${params.accessorydir} accessoryFiles
  
  cp accessoryFiles/scripts/genericFunctions.R .

  cp ${allele_AorCtypes} atype_ctype.txt
  
  R --no-save <${params.accessorydir}/scripts/drawAllevaFigure2.R
  """
  }

process drawAssociationsFigure{
  publishDir "${params.outdir}/figures",    mode: 'copy', overwrite: true, pattern: '*png'
  publishDir "${params.outdir}/figures",    mode: 'copy', overwrite: true, pattern: '*pdf'

  input:
  path(prdm9)
  path(cmh)
  path(type)
  path(clust)
  path(map)
  path(snpseq)

  output:
  path('Alleva_et_al_Fig*p*',     emit: mainFig)
  path('Alleva_et_al_Sup*p*',     emit: suppFigs)
  path('Alleva_et_al_Associ*tab', emit: cmhTable)

  """
  ln -s ${params.accessorydir} accessoryFiles

  for cmh in `ls *cmh`; do
    head \$cmh -n1 |perl -lane 'chomp; print \$_."\\tphenotype"' >header.tab
    phenotype=`echo \$cmh |perl -pi -e 's/prdm9AS\\.(.+?)\\.mod\\.cmh\$/\$1/' |perl -pi -e 's/like/\\-type/'`
    perl -lane 'print join("\\t",@F,'\$phenotype')' \$cmh >>table.tab
  done
  
  sort -k18,18 -k8rn,8rn table.tab |cut -f1-12,14-16,18 >table.sorted.tab
  
  cat header.tab table.sorted.tab >Alleva_et_al_Associtions.tab
  
  cp accessoryFiles/scripts/genericFunctions.R .
  cp accessoryFiles/scripts/drawAllevaAssociationPlots.R .
  
  R --no-save <drawAllevaAssociationPlots.R
  rm -f Rplots.pdf 
  """
  }
  
process drawRelatednessFigure{
  publishDir "${params.outdir}/figures",    mode: 'copy', overwrite: true, pattern: '*png'
  publishDir "${params.outdir}/figures",    mode: 'copy', overwrite: true, pattern: '*pdf'

  input:
  path(relatedness_table)
  path(allele_dets_A)
  path(allele_details)
  path(jeffreys_BS_alleles)
  path(pop_genotypes)

  output:
  path('Alleva_et_al_RelatednessPlot.p*',    emit: mainFig)
  path('Alleva_et_al_RelatednessPlot_S*.p*', emit: suppFig)

  """
  ln -s ${params.accessorydir} accessoryFiles
  
  ## Make table for Jeffreys data
  perl ${params.accessorydir}/scripts/makePRDM9variantsDF.pl 
  
  R --no-save <${params.accessorydir}/scripts/drawAllevaTemplateSwitchPlots.R
  
  rm -f Rplots.pdf
  """
  }

process drawGelQuantificationFigure{
  publishDir "${params.outdir}/figures",    mode: 'copy', overwrite: true, pattern: '*png'
  publishDir "${params.outdir}/figures",    mode: 'copy', overwrite: true, pattern: '*pdf'

  input:
  path(gtData)

  output:
  path('Allev*.p??',                  emit: mmFig)

  """
  ln -s ${params.accessorydir} accessoryFiles
  if [ -d "${params.outdir}" ]; then
    ln -s ${params.outdir}/final/dets/*details.txt .
    ln -s ${params.outdir}/fasta/rawmerged/*fa .
  else
    ln -s ${launchDir}/${params.outdir}/final/dets/*details.txt .
    ln -s ${launchDir}/${params.outdir}/fasta/rawmerged/*fa .
  fi

  ## GET ZFA sizes
  echo -e "id\\tpop\\tsize" >allZFAsizes.txt

  for x in *.prdm9_zf_array*details.txt; do
    name=`echo \${x/.prdm9_zf_array*/} | perl -pi -e 's/_/","/' 2>/dev/null`
    cat \$x |perl -lane 'if (\$_ =~ /ZF\\sCOUNT:\\s+(\\d+)/){print join("\\t","'\$name'",\$1)}' >>allZFAsizes.txt
  done

  ## Get Read lengths from fasta files
  echo -e "id\\tpop\\treadname\\tlength" >allReadLengths.txt

  for fa in *raw.fa; do
    name=`echo \${fa/.pb_ont*/} | perl -pi -e 's/_/","/' 2>/dev/null`
    perl -M"Bio::SeqIO" -se '\$fa = Bio::SeqIO->new(-file => \$fa);
                             while (\$s = \$fa->next_seq){
                              print join("\\t","'\$name'",\$s->display_id,length(\$s->seq)."\\n")}' -- -fa="\$fa" >>allReadLengths.txt
  done

  cp accessoryFiles/otherdata/*Gel*jpg .

  R --no-save <${params.accessorydir}/scripts/quantifyGels.R
  rm -f Rplots.pdf
  """
  }

process mergePublishedAndFoundAlleles{
  publishDir "${params.outdir}/annotation", mode: 'copy', overwrite: true, pattern: 'PrZFA*txt'

  input:
  path(dets)
  path(pubZFs)
  path(pubAlleles)

  output:
  path('PrZFA_alleles.txt',             emit: allele)
  path('PrZFA_ZFs.txt',                 emit: zf)
  path('distinct*txt',                  emit: names)

  """
  ln -s ${params.accessorydir} accessoryFiles

  perl -lane 'my @zfA = split /(.{84})/g,\$F[10];
              my @zfB = split /(.{84})/g,\$F[11];
              my @zfNmA = split /(.{2})/g,\$F[6];
              my @zfNmB = split /(.{2})/g,\$F[7];
              for my \$i(0..\$#zfA){
                if (\$zfNmA[\$i] && \$zfNmA[\$i] =~ /^\\!/){
                  print join("\\t",\$zfNmA[\$i],"NA","This_study",\$zfA[\$i]);
                }

                if (\$zfNmB[\$i] && \$zfNmB[\$i] =~ /^\\!/){
                  print join("\\t",\$zfNmB[\$i],"NA","This_study",\$zfB[\$i])
                }
              }' ${dets} |sort |grep -vP '^(NA|UU|NU)' |grep -vP '^\\s+\$' |uniq >PrZFA_ZFs.foundInPops.txt

  cat ${dets} |sort -V |grep -v allele |perl -lane 'print join("\\t",\$F[4],\$F[6]); print join("\\t",\$F[5],\$F[7])' |
               grep -vP '^(Unk|NA)' |sort |uniq  >PrZFA_alleles.foundInPops.txt

  sort -V ${pubZFs}     PrZFA_ZFs.foundInPops.txt     |uniq >PrZFA_ZFs.txt
  sort -V ${pubAlleles} PrZFA_alleles.foundInPops.txt |uniq >PrZFA_alleles.tmp

  ## Add logical fields for published alleles
  cut -f1 ${pubAlleles} >pubAlleles.IDs
  perl -lane 'BEGIN{open IN, "pubAlleles.IDs"; while (<IN>){chomp; \$ok{\$_}++};close IN}; 
              chomp; 
              print join("\\t",@F,\$ok{\$F[0]}?"TRUE":"FALSE")' PrZFA_alleles.tmp >PrZFA_alleles.vPub.tmp

  ## Add logical fields for pop alleles
  cut -f1 PrZFA_alleles.foundInPops.txt >PrZFA_alleles.foundInPops.IDs
  perl -lane 'BEGIN{open IN, "PrZFA_alleles.foundInPops.IDs"; while (<IN>){chomp; \$ok{\$_}++};close IN}; 
              chomp; 
              print join("\\t",@F,\$ok{\$F[0]}?"TRUE":"FALSE")' PrZFA_alleles.vPub.tmp >PrZFA_alleles.txt

  cut -f1 PrZFA_alleles.txt |uniq >distinct_PrZFA_allele_names.txt  

  """
  }

process inferRelatednessOfAlleles{
  input:
  path(alleles_to_check)
  path(PrZFA_alleles)
  //path(justForADelay)

  output:
  path('PrZFA_relatedness*tab', emit: pub)

  script:
  def randID = new Random().with {(1..30).collect {(('a'..'z')).join()[ nextInt((('a'..'z')).join().length())]}.join()}
  """
  ln -s ${params.accessorydir} accessoryFiles

  cp ${params.accessorydir}/scripts/checkForRecombinants.pl .
  cp ${params.accessorydir}/scripts/checkPrdm9BiParentalRecombinants.pl .
  
  perl checkForRecombinants.pl --allele ${alleles_to_check} >PrZFA_relatedness_graph.${randID}.tab
  """
  }

process getLinkedAlleles{
  publishDir "${params.outdir}/linkage", mode: 'copy', overwrite: true

  input:
  path(pr_haps)
  path(ACtypes)
  
  output:
  path('PRDM9*.bed'        , emit: prdm9)
  path('*.mod.*cmh'        , emit: cmh)
  path('*eno*txt'          , emit: type)
  path('*ust*txt'          , emit: clust)
  path('*eepers.sel*map'   , emit: map)
  path('*.SNPsequences.tab', emit: snpseq)

  script:
  """
  ln -s ${params.accessorydir} accessoryFiles
  
  ## Get position of PRDM9 gene & exons
  #wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.basic.annotation.gff3.gz
  #zcat gencode.v38.basic.annotation.gff3.gz |grep PRDM9 |grep -P '\\stranscript\\s' |cut -f1,4,5 |mergeBed -i -      >PRDM9.transcript.bed
  #zcat gencode.v38.basic.annotation.gff3.gz |grep PRDM9 |grep -P '\\sexon\\s'       |cut -f1,4,5 |sort -k1,1 -k2n,2n >PRDM9.exons.bed

  cp accessoryFiles/PRDM9locus/*bed .
  
  ## Hard code hg38 PRDM9 ZF locus coordinates
  zfstart=23526671
  zfend=23527765

  echo -e "chr5\\t\$zfstart\\t\$zfend" >PRDM9.ZFDomains.bed
  f=\$((zfstart - 20000000))
  t=\$((zfstart + 20000000))
  ## Change to whole chromosome
  #f=1
  #t=181538259
  
  ## Download 1KG SNPs ... tabix slicing not working for some reason. 
  wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.chr5.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz
  wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.chr5.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz.tbi
  #gsutil cp gs://genomics-public-data/1000-genomes/vcf/ALL.chr5.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.vcf .
  #gsutil cp gs://genomics-public-data/1000-genomes/vcf/ALL.chr5.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.vcf.tbi .
  
  ## Re-label as VCFv4.2 : vcftools not compatible with v4.3 !!
  #  NOTE: this just tricks vcftools - format doesn't change, but for this purpose, that's OK!
  tabix -h ALL.chr5.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz 5:\$f-\$t |perl -pi -e 's/VCFv4.3/VCFv4.2/' >reg.vcf
  
  ## Get all individuals with a called diploid PRDM9 genotype
  ## exclude individuals with an unknown allele
  ## exclude YRI "children" - both parents also in cohort 
  grep -v ^id ${pr_haps} |grep -vP '(Unk|NoZFA)' |grep -v TRUE\$ >pr_haps_to_use.txt
  cat pr_haps_to_use.txt |cut -f1 |sort |uniq >individuals.txt
  
  ## Keep genotypes only for select individuals (with diploid Prdm9 genotype)
  vcftools --keep individuals.txt --gzvcf reg.vcf \
          --from-bp \$f --chr 5 --to-bp \$t \
          --mac 1 --max-mac 10000 \
          --recode --remove-indels \
          --maf 0.02 \
          --max-missing 1 \
          --out prdm9AS_keepers

  ## Get final list of people to use and count
  head -n 2000 prdm9AS_keepers.recode.vcf |grep ^# |tail -n1 |\
                        perl -pi -e 's/^.+INFO\\s+FORMAT\\s+//' 2>e.e |\
                        perl -lane 'print join("\\n",@F)' >individuals_for_plink.txt
                        
  ## Add rs IDS to 1KG vcf file
  #tabix -h https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.38.gz NC_000005.10:\$f-\$t |\
  #      perl -pi -e 's/NC_000005.10/5/g' >regionRS.vcf
  
  tabix -h https://ftp.ncbi.nih.gov/snp/archive/b154/VCF/GCF_000001405.38.gz  NC_000005.10:\$f-\$t |\
        perl -pi -e 's/NC_000005.10/5/g' >regionRS.vcf
  
  ncols=`head -n 2000 prdm9AS_keepers.recode.vcf |grep ^# |tail -n1 |perl -lane 'print \$#F+1'`

  intersectBed -a prdm9AS_keepers.recode.vcf \
               -b regionRS.vcf -wao |\
               perl -lane 'if (\$_ =~ /^.+\\s(rs\\d+)/){\$F[2] = \$1}; \$out=join("\\t",@F); print \$out' |\
               cut -f1-\$ncols |\
               sort -k1,1 -k2n,2n |\
               uniq >prdm9AS_keepers.tmpvcf

  grep ^# prdm9AS_keepers.recode.vcf >head.vcf
  cat head.vcf prdm9AS_keepers.tmpvcf >prdm9AS_keepers.vcf
  
  bgzip prdm9AS_keepers.vcf
  tabix -p vcf prdm9AS_keepers.vcf.gz 

  ## Make plink files for "final" VCF
  vcftools --gzvcf prdm9AS_keepers.vcf.gz \
           --plink --out prdm9AS_keepers

  ## Get phenotype-individual table
  ## Add A-type/C-type to table
  grep -v ^# PrZFA_alleles.details.txt |sort >ACsorted.txt
  
  grep -f individuals_for_plink.txt ${pr_haps} >haps_plink.1.tmp
  
  sort -k6,6 haps_plink.1.tmp >haps_plink.2.tmp

  join -1 6 -2 1 haps_plink.2.tmp ACsorted.txt |perl -pi -e 's/\\s+(\\S)/\\t\$1/g' |perl -pi -e 's/\\s+(\\S)/\\t\$1/g' |cut -f1-7,25 |sort -k6,6 >haps_plink.3.tmp
  join -1 6 -2 1 haps_plink.3.tmp ACsorted.txt |perl -pi -e 's/\\s+(\\S)/\\t\$1/g' |perl -pi -e 's/\\s+(\\S)/\\t\$1/g' |cut -f1-8,18             >prdm9_haplotypes.tmp

  perl -lane 'print join("\\t",@F[2..5],@F[0..1],\$F[7]."like",\$F[8]."like")' prdm9_haplotypes.tmp >prdm9_haps_for_plink.tab

  ## Build phenotypes, clusters (populations) and genotypes files for plink
  cat prdm9_haps_for_plink.tab |\
      perl -lane 'BEGIN{print join("\\t","FID","IID","A","B","C","L14","Alike","Clike")}; \
                        \$A = (\$F[4] eq "A" || \$F[5] eq "A")?2:1; 
                        \$B = (\$F[4] eq "B" || \$F[5] eq "B")?2:1; 
                        \$C = (\$F[4] eq "C" || \$F[5] eq "C")?2:1; 
                        \$L14 = (\$F[4] eq "L14" || \$F[5] eq "L14")?2:1; 
                        \$Alike = (\$F[6] eq "Alike" || \$F[7] eq "Alike")?2:1; 
                        \$Clike = (\$F[6] eq "Clike" || \$F[7] eq "Clike")?2:1; 
                        print join("\\t",\$F[0],\$F[0],\$AA,\$A,\$B,\$C,\$L14,\$Alike,\$Clike)' >phenotypes.txt
                        
  grep -v ^id prdm9_haps_for_plink.tab |awk '{print \$1 "\\t" \$1 "\\t" \$2}' >clusters.txt
  grep -v ^id prdm9_haps_for_plink.tab |awk '{print \$1 "\\t" \$1 "\\t" \$3}' >genotypes.txt

  ## Run association tests
  plink --file prdm9AS_keepers --assoc --pheno phenotypes.txt --all-pheno --allow-no-sex --within clusters.txt --mh

  ## Generate output files for individual phenotypes
  for phen in "A" "B" "C" "L14" "Alike" "Clike"; do
    
    ## Label SNPs that coincide with PRDM9 transcript, ZFDomains, Exons
    perl -lane 'print join("\\t","chr5",\$F[2]-1,\$F[2],\$F[1]) unless (\$_ =~ /CHR/)' plink.\$phen.cmh |\
                intersectBed -a - -b PRDM9.transcript.bed -c |\
                intersectBed -a - -b PRDM9.ZFDomains.bed  -c |\
                intersectBed -a - -b PRDM9.exons.bed  -c |\
                perl -lane 'BEGIN{print join("\\t","rsagain","transcript","zfdomain","exon")}; 
                                  print join("\\t",@F[3..6])' >ol.tab

    paste plink.\$phen.cmh ol.tab >plinkMOD.\$phen.cmh

    grep -v CHR plinkMOD.\$phen.cmh |\
                perl -lane 'if (\$F[7] !~ /^[0123456789\\.e\\-]+\$/ || \$F[7] eq "0" || \$F[14] > 0){
                              \$F[7] = -1;
                            }else{
                              \$F[7]=-1*log(\$F[7])/log(10);
                            }
                            print join("\\t",@F)' |\
                sort -k8rn,8rn >associated_SNPs.\$phen.tab

    head -n 10       associated_SNPs.\$phen.tab  >selectSNPs.\$phen.tmp
    grep rs6889665   associated_SNPs.\$phen.tab >>selectSNPs.\$phen.tmp
    grep rs1603084   associated_SNPs.\$phen.tab >>selectSNPs.\$phen.tmp
    grep rs10057021  associated_SNPs.\$phen.tab >>selectSNPs.\$phen.tmp
    grep rs141586808 associated_SNPs.\$phen.tab >>selectSNPs.\$phen.tmp
    grep rs139754603 associated_SNPs.\$phen.tab >>selectSNPs.\$phen.tmp

    sort -k1,1 selectSNPs.\$phen.tmp |uniq >selectSNPs.\$phen.tab
    
    cat selectSNPs.\$phen.tab |\
                perl -pi -e 's/\\s+(\\S)/\\t\$1/g' |\
                cut -f2 >selectSNPs.\$phen.list

    grep -vf selectSNPs.\$phen.list plinkMOD.\$phen.cmh |perl -lane 'chomp; print \$_."\\t".(\$line++?"0":"selected")' >others.txt
    grep -f  selectSNPs.\$phen.list plinkMOD.\$phen.cmh |perl -lane 'chomp; print \$_."\\t1"'                          >select.txt
    sort -k3n,3n others.txt select.txt |perl -pi -e 's/^\\s+//' |perl -pi -e 's/\\s+(\\S+)/\\t\$1/g'         >prdm9AS.\$phen.mod.cmh

    rm others.txt select.txt
    
    ## Make genotypes file for select SNPs
    perl -lane 'print join("\\t",\$F[0],\$F[2]-1,\$F[2])' selectSNPs.\$phen.tab >selectSNPs.\$phen.bed
    tabix -h prdm9AS_keepers.vcf.gz -R selectSNPs.\$phen.bed >selectSNPs.\$phen.vcf

    vcftools --vcf selectSNPs.\$phen.vcf \
             --plink \
             --out prdm9AS_keepers.select\$phen
    
    ## Change order in ped file
    headerTxt="FID"
    for x in `grep -vP '^#' selectSNPs.\$phen.vcf |cut -f3`; do 
      headerTxt=\$headerTxt"\\t"\$x"_A\\t"\$x"_B"
    done
    
    echo -e \$headerTxt >prdm9AS_keepers.select\$phen.SNPsequences.tab
    cut -f1,7-1000 prdm9AS_keepers.select\$phen.ped >>prdm9AS_keepers.select\$phen.SNPsequences.tab
    
    
  done


  """
  }
  
process dnaToPeptide{
  publishDir "${params.outdir}/annotation", mode: 'copy', overwrite: true, pattern: '*txt'

  input:
  path(PrZFA_alleles)
  path(PrZFA_ZFs)
  path(allBloodAndSperm_alleles)
  
  output:
  path('PrZFA_alleles.details.txt', emit: alleles)
  path('PrZFA_ZFs.details.txt'    , emit: zfs)

  script:
  """
  ln -s ${params.accessorydir} accessoryFiles

  ## get allele NT sequences
  perl -lane 'BEGIN{open IN, "${PrZFA_ZFs}";
              while (<IN>){
                  my @X = split(/\\t/,\$_);
                  chomp \$X[3];
                  \$aa{\$X[0]} = \$X[3];
                };
                close IN;
              };
              my \$seq;
              chomp \$F[3];
              \$F[0] =~ s/^\\s*\\>(\\S+)/\$1/;
              while (\$F[1] =~ s/^(\\S\\S)//){\$seq .= \$aa{\$1}};
              print \$F[0]."\\t".\$seq' ${PrZFA_alleles} |sort >PrZFA_alleles.NT.txt

  awk '{print ">"\$1"\\n"\$2}' PrZFA_alleles.NT.txt |fold -w 28 >PrZFA_alleles.NT.fa

  ## ZFs to AA sequences
  perl -lane 'open RAW, ">", "x.raw";
              print RAW \$F[3];
              close RAW;
              \$aa=`transeq -sequence x.raw -outseq /dev/stdout 2>/dev/null |grep -vP "^\\>"`;
              chomp \$aa;
              print \$F[0]."\\t".\$aa' ${PrZFA_ZFs} |sort >PrZFA_ZFs.AA.txt

  awk '{print ">"\$1"\\n"\$2}' PrZFA_ZFs.AA.txt |fold -w 28 >PrZFA_ZFs.AA.fa

  ## Alleles to ZF sequences
  perl -lane 'BEGIN{open IN, "PrZFA_ZFs.AA.txt";
              while (<IN>){
                  my @X = split(/\\t/,\$_);
                  chomp \$X[1];
                  \$aa{\$X[0]} = \$X[1];
                };
                close IN;
              };
              my \$AAseq;
              chomp \$F[3];
              \$F[0] =~ s/^\\s*\\>(\\S+)/\$1/;
              while (\$F[1] =~ s/^(\\S\\S)//){\$AAseq .= \$aa{\$1}};
              print \$F[0]."\\t".\$AAseq' ${PrZFA_alleles} |sort >PrZFA_alleles.AA.txt

  awk '{print ">"\$1"\\n"\$2}' PrZFA_alleles.AA.txt |fold -w 28 >PrZFA_alleles.AA.fa

  ## ZFs to contact residues
  perl -lane 'open RAW, ">", "x.raw";
              print RAW \$F[3];
              close RAW;
              \$aa=`transeq -sequence x.raw -outseq /dev/stdout 2>/dev/null |grep -vP "^\\>"`;
              chomp \$aa;
              @AA = split(//,\$aa);
              print \$F[0]."\\t".join("",\$AA[6],\$AA[7],\$AA[9],\$AA[12])' ${PrZFA_ZFs} |sort >PrZFA_ZFs.AAcontactresidues.txt

  awk '{print ">"\$1"\\n"\$2}' PrZFA_ZFs.AAcontactresidues.txt |fold -w 28 >PrZFA_ZFs.AAcontactresidues.fa

  ## Alleles to contact residues
  perl -lane 'BEGIN{open IN, "PrZFA_ZFs.AAcontactresidues.txt";
              while (<IN>){
                  my @X = split(/\\t/,\$_);
                  chomp \$X[1];
                  \$aa{\$X[0]} = \$X[1];
                };
                close IN;
              };
              my \$AAseq;
              chomp \$F[3];
              \$F[0] =~ s/^\\s*\\>(\\S+)/\$1/;
              while (\$F[1] =~ s/^(\\S\\S)//){\$AAseq .= \$aa{\$1}};
              print \$F[0]."\\t".\$AAseq' ${PrZFA_alleles} |sort >PrZFA_alleles.AAcontactresidues.txt

  awk '{print ">"\$1"\\n"\$2}' PrZFA_alleles.AAcontactresidues.txt |fold -w 28 >PrZFA_alleles.AAcontactresidues.fa

  #############################################################################
  ## Now collect these into one big TSV for alleles
  ## Do ZF info first
  echo -e "#code           The two-character ZF code used in this publication"     >PrZFA_ZFs.details.txt
  echo -e "#old_code       The one-character ZF code used in past publication(s)" >>PrZFA_ZFs.details.txt
  echo -e "#source         The original source that described this ZF"            >>PrZFA_ZFs.details.txt
  echo -e "#dna_sequence   The DNA sequence of this ZF"                           >>PrZFA_ZFs.details.txt
  echo -e "#aa_sequence    The amino acid sequence of this ZF"                    >>PrZFA_ZFs.details.txt
  echo -e "dna_contact_aas The amino acids that contact DNA (-1,2,3,6 positions; see https://doi.org/10.1074/jbc.M117.805754)" >>PrZFA_ZFs.details.txt

  echo -e "code\told_code\tsource\tdna_sequence\taa_sequence\tdna_contact_aas"     >>PrZFA_ZFs.details.txt

  sort -k1,1 PrZFA_ZFs.txt >sortedZFs.txt
  join sortedZFs.txt PrZFA_ZFs.AA.txt                |perl -pi -e 's/\\s+(\\S)/\\t\$1/g' |sort >zf1.txt
  join zf1.txt       PrZFA_ZFs.AAcontactresidues.txt |perl -pi -e 's/\\s+(\\S)/\\t\$1/g' |sort >>PrZFA_ZFs.details.txt

  #############################################################################
  ## Check A-type / C-type alleles
  #############################################################################

  grep -P '^A\\s' PrZFA_alleles.AAcontactresidues.txt |perl -lane '\$_ =~ /^.+(.{24})\$/; print ">A-BD\\n\$1"' >A-BD.fa
  grep -P '^C\\s' PrZFA_alleles.AAcontactresidues.txt |perl -lane '\$_ =~ /^.+(.{28})\$/; print ">C-BD\\n\$1"' >C-BD.fa

  makeblastdb -dbtype prot -in PrZFA_alleles.AAcontactresidues.fa -out PRZFA

  blastp -query A-BD.fa -db PRZFA -outfmt 7 -gapopen 32767 -gapextend 32767 -max_hsps 1 -num_alignments 10000 >A-BD_hits.blastp.txt
  blastp -query C-BD.fa -db PRZFA -outfmt 7 -gapopen 32767 -gapextend 32767 -max_hsps 1 -num_alignments 10000 >C-BD_hits.blastp.txt

  cat A-BD_hits.blastp.txt |grep -v ^# |perl -lane 'print \$F[1]."\\t".(((24-\$F[3]))+\$F[4])' |sort -k1,1 >ABDhits.txt
  cat C-BD_hits.blastp.txt |grep -v ^# |perl -lane 'print \$F[1]."\\t".(((28-\$F[3]))+\$F[4])' |sort -k1,1 >CBDhits.txt

  echo -e "allele\\tdA\\tdC"                  >AC_BDhits.txt
  paste ABDhits.txt CBDhits.txt |cut -f1,2,4 >>AC_BDhits.txt

  paste ABDhits.txt CBDhits.txt |cut -f1,2,4 |perl -lane 'BEGIN{print "allele\ttype"};
                                              print join("\\t",\$F[0], \$F[1], \$F[2], (\$F[1]<\$F[2]?"A":"C"))' |sort >PrZFA_alleles.ACtype.txt

  echo -e "#ID               PRDM9 allele long ID"                                     >PrZFA_alleles.details.txt
  echo -e "#short_ID         PRDM9 allele short ID (for alleles from Jeffreys 2013)"  >>PrZFA_alleles.details.txt
  echo -e "#published_allele TRUE if allele is from a previously published work"      >>PrZFA_alleles.details.txt
  echo -e "#in_pop           TRUE if allele is found in human populations"            >>PrZFA_alleles.details.txt
  echo -e "#                 FALSE if only found in sperm / blood from Jeffreys 2013" >>PrZFA_alleles.details.txt
  echo -e "#dna_sequence     The DNA sequence of the ZF array for this allele"        >>PrZFA_alleles.details.txt
  echo -e "#aa_sequence      The amino acid sequence of the ZF array for this allele" >>PrZFA_alleles.details.txt
  echo -e "#dna_contact_aas  The amino acids that contact DNA (-1,2,3,6 positions; see https://doi.org/10.1074/jbc.M117.805754)"   >>PrZFA_alleles.details.txt
  echo -e "#ABDdist          BLAST distance to PRDM9A binding site"                   >>PrZFA_alleles.details.txt
  echo -e "#CBDdist          BLAST distance to PRDM9C binding site"                   >>PrZFA_alleles.details.txt
  echo -e "#ACtype           Is this allele more A-type or C-type"                    >>PrZFA_alleles.details.txt

  echo -e "ID\\tshort_ID\\tzf_code\\tpublished_allele\\tin_pop\\tdna_sequence\\taa_sequence\\tdna_contact_aas\\tABDdist\\tCBDdist\\tACtype" >>PrZFA_alleles.details.txt

  ## Sort & add shortID column
  sort -k1,1 PrZFA_alleles.txt | perl -lane '\$sID = \$F[0];
                                             \$sID =~ s/^(\\S+v):\\S+?:(\\d+):.+\$/\$1:\$2/;
                                             \$F[0] = \$F[0]."\\t".\$sID;
                                             print join("\\t",@F)' >sortedalleles.txt

  join sortedalleles.txt PrZFA_alleles.NT.txt                |perl -pi -e 's/\\s+(\\S)/\\t\$1/g' |sort >s1.txt
  join s1.txt            PrZFA_alleles.AA.txt                |perl -pi -e 's/\\s+(\\S)/\\t\$1/g' |sort >s2.txt
  join s2.txt            PrZFA_alleles.AAcontactresidues.txt |perl -pi -e 's/\\s+(\\S)/\\t\$1/g' |sort >s3.txt
  join s3.txt            PrZFA_alleles.ACtype.txt            |perl -pi -e 's/\\s+(\\S)/\\t\$1/g' |sort >>PrZFA_alleles.details.txt

  """
  }

process makeACtypesPlot{
  publishDir "${params.outdir}/figures",    mode: 'copy', overwrite: true, pattern: '*png'
  publishDir "${params.outdir}/figures",    mode: 'copy', overwrite: true, pattern: '*pdf'

  input:
  path(allele_dets)
  path(zf_dets)
  path(allBloodAndSperm_alleles)
  
  output:
  path('*.p??', emit: img)

  script:
  """
  ln -s ${params.accessorydir} accessoryFiles

  #perl -lane 'print ">\$F[1]\\n\$F[6]" if (\$F[4] eq "TRUE")' ${allele_dets} |fold -w 28 >pubAlleles.AA.fa
  #pwm_predict -m SVMp pubAlleles.AA.fa
  
  cp ${allele_dets} atype_ctype.txt
  
  cp accessoryFiles/otherdata/pubAlleles.AA.pwm .
  cp accessoryFiles/otherdata/PrZFA_Alleles.AA.pwm .
  
  cp accessoryFiles/img/motifAlignmentManual.png .
  
  R --no-save <accessoryFiles/scripts/drawACtypePlot.R
  """  
  }

process processHumanHotspots {

  publishDir "${params.outdir}/humanHS/table", mode: 'copy', overwrite: true, pattern: '*tab'
  publishDir "${params.outdir}/humanHS/fig",   mode: 'copy', overwrite: true, pattern: '*png'
  publishDir "${params.outdir}/humanHS/fasta",   mode: 'copy', overwrite: true, pattern: '*fa'
  //publishDir "output/humanHS/fig",   mode: 'copy', overwrite: true, pattern: '*pdf'

  input:
  path(hotspots)

  output:
  path('hotspotsData.tab',  emit: table)
  path('*png',              emit: img)
  path('*chipmunk.fa',      emit: chipmunkFA)
  path('*meme.fa',          emit: memeFA)
  path('*seq.fa',           emit: fasta)
  path('*bigwin.fa',        emit: wideFA)

  """
  ln -s ${params.accessorydir} accessoryFiles

  ## Merge all HS
  for hs in *.peaks.bedgraph; do
    new=\${hs/.peaks.bedgraph/.400bp.bedgraph}
    slopBed -g ${genomeIDX} -i \$hs -l -0.5 -r -0.5 -pct |slopBed -g ${genomeIDX} -i - -l 200 -r 200 |sort -k1,1 -k2n,2n >\$new

    hottest=\${hs/.peaks.bedgraph/.hottest.bed}
    cat \$hs |sort -k4rn,4rn |head -n 100 |sort -k1,1 -k2n,2n |cut -f1-3 >\$hottest
  done

  sort -k1,1 -k2n,2n -k3n,3n *400bp.bedgraph |mergeBed -i - |cut -f1-3 >all.bed

  for hs in *400bp.bedgraph; do
    sort -k1,1 -k2n,2n -k3n,3n \$hs -o \$hs

    name=\${hs/.400bp.bedgraph/}
    echo \$name >\$name".OL"

    mapBed -a all.bed -b \$hs -c 4 -o sum |perl -pi -e 's/\\s+\\./\\t0/' |cut -f4 >>\$name".OL"
  done

  echo -e "cs\tfrom\tto" >hotspotsData.tmp
  cat all.bed >>hotspotsData.tmp

  paste hotspotsData.tmp *OL >hotspotsData.tab

  cp accessoryFiles/scripts/genericFunctions.R .
  cp accessoryFiles/scripts/plot_AA_v_AA_v_AN_comparison.R .
  cp accessoryFiles/scripts/getDifferentialHS.R .

  R --vanilla <plot_AA_v_AA_v_AN_comparison.R
  R --vanilla <getDifferentialHS.R

  sort -k4rn L4_biased_HS.bedgraph |head -n ${params.npeaks} |sort -k1,1 -k2n,2n |cut -f1-3 >L4_biased_Auto.top.bed
  sort -k4rn B_biased_HS.bedgraph  |head -n ${params.npeaks} |sort -k1,1 -k2n,2n |cut -f1-3 >B_biased_Auto.top.bed
  sort -k4rn N_biased_HS.bedgraph  |head -n ${params.npeaks} |sort -k1,1 -k2n,2n |cut -f1-3 >N_biased_Auto.top.bed

  ## Get differential hotspots
  #grep N-up MAdetails_AA_v_AN.tab   |cut -f1-3,10 |sort -k1,1 -k2n,2n >AN_comp_Nup.bedgraph
  #grep N-down MAdetails_AA_v_AN.tab |cut -f1-3,8  |sort -k1,1 -k2n,2n >AN_comp_Ndown.bedgraph
  #grep Stable MAdetails_AA_v_AN.tab |cut -f1-3,8  |sort -k1,1 -k2n,2n >AN_comp_stable.bedgraph
  #grep N-up MAdetails_AA_v_AN.tab   |cut -f1-3    |sort -k1,1 -k2n,2n >AN_comp_Nup.bed
  #grep N-down MAdetails_AA_v_AN.tab |cut -f1-3    |sort -k1,1 -k2n,2n >AN_comp_Ndown.bed
  #grep Stable MAdetails_AA_v_AN.tab |cut -f1-3    |sort -k1,1 -k2n,2n >AN_comp_stable.bed

  #grep N-up MAdetails_AA_v_AN.tab   |cut -f1-3,10 |sort -k4rn |head -n ${params.npeaks} |sort -k1,1 -k2n,2n |cut -f1-4 >AN_comp_Nup_topALL.bedgraph
  #grep N-down MAdetails_AA_v_AN.tab |cut -f1-3,8  |sort -k4rn |head -n ${params.npeaks} |sort -k1,1 -k2n,2n |cut -f1-4 >AN_comp_Ndown_topALL.bedgraph
  #grep Stable MAdetails_AA_v_AN.tab |cut -f1-3,8  |sort -k4rn |head -n ${params.npeaks} |sort -k1,1 -k2n,2n |cut -f1-4 >AN_comp_stable_topALL.bedgraph
  #grep N-up MAdetails_AA_v_AN.tab   |cut -f1-3,10 |sort -k4rn |head -n ${params.npeaks} |sort -k1,1 -k2n,2n |cut -f1-3 >AN_comp_Nup_topALL.bed
  #grep N-down MAdetails_AA_v_AN.tab |cut -f1-3,8  |sort -k4rn |head -n ${params.npeaks} |sort -k1,1 -k2n,2n |cut -f1-3 >AN_comp_Ndown_topALL.bed
  #grep Stable MAdetails_AA_v_AN.tab |cut -f1-3,8  |sort -k4rn |head -n ${params.npeaks} |sort -k1,1 -k2n,2n |cut -f1-3 >AN_comp_stable_topALL.bed

  ##Get top N hotspots for motif finding
  npeaks=${params.npeaks}
  bigwinSize=${params.bigwinSize}
  centerSz=\$((${params.centerSz}/2))
  midpoint=\$((centerSz/2+1))

  sort -k4n,4n AA1.400bp.bedgraph |head -n \$npeaks |cut -f1-3 |sort -k1,1 -k2n,2n >A_top.bed

  sort -k1,1 -k2n,2n AA*400bp.bedgraph |cut -f1-3 >allAA.bed
  intersectBed -a AA1.400bp.bedgraph -b allAA.bed -u |sort -k1,1 -k2n,2n >A_all.bedgraph
  intersectBed -a AB1.400bp.bedgraph -b allAA.bed -v |sort -k1,1 -k2n,2n >B_all.bedgraph
  intersectBed -a AN.400bp.bedgraph  -b allAA.bed -v |sort -k1,1 -k2n,2n >N_all.bedgraph

  intersectBed -a AC.400bp.bedgraph   -b allAA.bed      -v |sort -k1,1 -k2n,2n >C_all.bedgraph
  intersectBed -a CL4.400bp.bedgraph  -b C_all.bedgraph -v |sort -k1,1 -k2n,2n >L4_all.bedgraph

  sort -k4rn,4rn A_all.bedgraph |head -n \$npeaks |cut -f1-3 |sort -k1,1 -k2n,2n >A_topALL.bed
  sort -k4rn,4rn N_all.bedgraph |head -n \$npeaks |cut -f1-3 |sort -k1,1 -k2n,2n >N_topALL.bed
  sort -k4rn,4rn C_all.bedgraph |head -n \$npeaks |cut -f1-3 |sort -k1,1 -k2n,2n >C_topALL.bed
  sort -k4rn,4rn B_all.bedgraph |head -n \$npeaks |cut -f1-3 |sort -k1,1 -k2n,2n >B_topALL.bed

  intersectBed -a CL4.400bp.bedgraph  -b C_all.bedgraph -v |sort -k4rn,4rn |head -n \$npeaks |cut -f1-3 |sort -k1,1 -k2n,2n >L4_topALL.bed

  for bed in *topALL.bed; do
    type=\${bed/_topALL.bed/}

    if [ `grep -P 'chr[XY]' \$bed  |wc -l ` -gt 1 ]; then
      grep -P 'chr[XY]'  \$bed >\$type"_XY.top.bed"
    fi

    if [ `grep -vP 'chr[XY]' \$bed  |wc -l ` -gt 1 ]; then
      grep -vP 'chr[XY]'  \$bed >\$type"_Auto.top.bed"
    fi

  done

  for bed in *Auto.top.bed; do
    fa=\${bed/bed/seq.fa}
    chipmunk_fa=\${bed/bed/chipmunk.fa}
    meme_fa=\${bed/bed/meme.fa}
    slopBed -i \$bed -l -0.5 -r -0.5 -g ${genomeIDX} -pct |slopBed -i - -g ${genomeIDX} -l \$centerSz -r \$centerSz >centers.bed
    bedtools getfasta -fi ${genomeFA} -bed centers.bed -fo \$fa

    if [ "${params.RM}" == "false" ]; then
      cat \$fa | tr [a-z] [A-Z] >uc.fa
      mv uc.fa \$fa
    fi

    cat \$fa |perl -lane 'if (\$_ =~ /^\\>/){print \$_}else{\$_ =~ s/[gatc]/N/g; print \$_}' |fold -w 60 >\$meme_fa
    cat \$fa |perl -lane 'if (\$_ =~ /^\\>/){print "\\>'\$midpoint'"}else{\$_ =~ s/[gatc]/N/g; print \$_}' |fold -w 60 >\$chipmunk_fa
  done

  for bg in *_all.bedgraph; do
    ## Make wide windows for motif finding later
    faBW=\${bg/_all.bedgraph/.bigwin.fa}
    tmpBG=\${bg/_all.bedgraph/_tmp1.bedgraph}
    renBG=\${bg/_all.bedgraph/_tmp2.bedgraph}
    tmpFA=\${bg/_all.bedgraph/_tmp.fa}
    tmpNm=\${bg/_all.bedgraph/}

    bgname=\${bg/_all.bedgraph}

    slopBed -i \$bg -l -0.5 -r -0.5  -g ${genomeIDX} -pct |slopBed -i - -g ${genomeIDX} -l \$bigwinSize -r \$bigwinSize |shuf |head -n 5000 |sort -k1,1 -k2n,2n >\$tmpBG
    perl -lane 'print join("\t",@F[0..2],join(":",@F,"'\$tmpNm'"))' \$tmpBG >\$renBG
    bedtools getfasta -name -fi ${genomeFA} -bed \$renBG -fo \$tmpFA
    ##cat \$tmpFA |perl -lane 'if (\$_ =~ /^\\>/){print \$_}else{\$_ =~ s/[gatc]/N/g; print \$_}' |fold -w 60 >\$faBW
    fold -w 60 \$tmpFA >\$faBW
  done

  """
  }

process getMotifsMEME {

  publishDir "${params.outdir}/humanHS/motifs/meme",        mode: 'copy', overwrite: true

  tag { fasta }

  input:
  path(fasta)
  path(wideFA)

  output:
  path("*motifdensity*", optional: true, emit: motifdist)
  path("*MEME*png",      optional: true, emit: png)
  path("*MEME*eps",      optional: true, emit: eps)
  path("meme_details*",  optional: true, emit: meme)

  script:
  def name     = fasta.name.replaceFirst(".meme.fa","")
  def maskedFA = fasta.name.replaceFirst(".meme.fa",".meme.masked.fa")
  def bigwin   = fasta.name.replaceFirst(".bigwin.fa","")

  """
  #for seed in {1..5}; do
  #cat ${fasta} |awk '/^>/ { if(i>0) printf("\\n"); i++; printf("%s\\t",\$0); next;} {printf("%s",\$0);} END { printf("\\n");}' >faOneLiner.txt
  #cat faOneLiner.txt |shuf  --random-source=<(cat faOneLiner.txt) |head -n 500 |awk '{printf("%s\\n%s\\n",\$1,\$2)}'

  dust ${fasta} >${maskedFA}
  meme ${maskedFA} -oc meme_out -dna -objfun ce -nmotifs 5 -revcomp -seed 42 -minw 12 -maxw 26 -minsites 100 -maxsize 1020000
  centrimo -oc centrimo_out ${bigwin} meme_out/meme.txt

  mfile=${name}".MEMEmotifdensity.txt";
  echo -e "hotspots\\ttype\\tmotifsrc\\tscore\\tpos\\tstrand\\tcs\\tfrom\\tto\\tstrength\\tallele" >\$mfile

  for motifLine in `cut -f2,3,5,6,7 centrimo_out/centrimo.tsv |grep MEME |perl -lane 'print join(":",\$F[0],\$F[1]) if (!\$cnt++ || \$F[2] < 1e-3)'`; do

    motif=`echo \$motifLine | cut -d":" -f1`
    motifName=`echo \$motifLine |cut -d":" -f2`

    echo \$motif" --- "\$motifName
    logos=\${motifName/MEME-/}

    cp "meme_out/logo"\$logos".png"    ${name}.MEME-\$motifName".FWD.png"
    cp "meme_out/logo_rc"\$logos".png" ${name}.MEME-\$motifName".REV.png"
    cp "meme_out/logo"\$logos".eps"    ${name}.MEME-\$motifName".FWD.eps"
    cp "meme_out/logo_rc"\$logos".eps" ${name}.MEME-\$motifName".REV.eps"

    for fa in *.bigwin.fa; do
      genotype=\${fa/.bigwin.fa};

      fimo --motif \$motif --oc fimo_out meme_out/meme.txt \$fa
      perl -lane '@X=split(/:/,\$F[0]); unless (\$_ =~ /version/){print join("\\t","'\$genotype'","meme","${name}",\$F[5],\$F[3],\$F[6],@X[0..3],"'\$genotype'")}' fimo_out/fimo.gff >>\$mfile
      mv fimo_out fimo_${name}_\$motifName"_"\$genotype"_out"
    done

  done

  mkdir meme_details_${name}
  cp -r meme_out     meme_details_${name}
  cp -r centrimo_out meme_details_${name}
  cp -r fimo_*_out   meme_details_${name}
  """
  }

process drawHotspotsFigure {

  publishDir "${params.outdir}/figures",          mode: 'copy', overwrite: true, pattern: '*png'
  publishDir "${params.outdir}/figures",          mode: 'copy', overwrite: true, pattern: '*pdf'
  publishDir "${params.outdir}/humanHS/motifPWM", mode: 'copy', overwrite: true, pattern: '*yaml'

  input:
  path(meme)
  path(hs)

  output:
  path("*png",                       optional: true, emit: png)
  path("*pdf",                       optional: true, emit: pdf)
  path("*predictedBindingSite*yaml", optional: true, emit: prdm9PredictedBS)
  path("*topAutosomal*yaml",         optional: true, emit: prdm9Motifs)

  script:
  """
  ln -s ${params.accessorydir} accessoryFiles

  cp accessoryFiles/otherdata/selectPRDM9.AA.fa .
  cp accessoryFiles/otherdata/selectPRDM9.AA.pwm .

  ## Predict PRDM9 binding
  ## pwm_predict isn't playing nice with docker container, so
  ## I just included the pwm file in accessorydata
  ## The command used is :
  #./pwm_predict -m SVMp selectPRDM9.AA.fa

  R --no-save <${params.accessorydir}/scripts/drawAllevaFigure4.R
  R --no-save <${params.accessorydir}/scripts/drawAllevaHotspotsFigure_Supp.R
  rm -f Rplots.pdf
  """
  }

workflow faToGenotype{
  take:
  name
  fa
  publishedZFs
  publishedAlleles

  main:
  just_ZF_fastas = extractZFsfromRawFA(name,
                                       fa.collate(10),
                                       publishedZFs,
                                       publishedAlleles)

  haps     = makeHaplotypeTable(name,
                                just_ZF_fastas.hap.collect(),
                                just_ZF_fastas.fa.collect(),
                                fa.collect())

  zflens   = genotypeVsZFlengthDistributions(name,
                                             haps.hap,
                                             just_ZF_fastas.dets.collect())

  emit:
  haps.hap

  }

workflow genotypeFromPacbio {
  take:
  pubZFs
  pubAlleles

  main:
  //raw_fasta = Channel.fromPath("${params.pbfa}/*.fa")
  pb_fa = Channel.fromPath("${params.pbfa}/*.fa")
  raw_fasta = merge_raw_fa(pb_fa.collect(),'pb')

  // This is required to convert the ArrayList output to a Channel !
  raw_fa     = raw_fasta.fa.flatten().map { sample -> file(sample) }

  genotypes = faToGenotype("pacbio", raw_fa,  pubZFs, pubAlleles)

  emit:
  genotypes
  }

workflow genotypeFromONT {
  take:
  pubZFs
  pubAlleles

  main:
  in_fa = Channel.fromPath("${params.ontfa}/*.fa")
  raw_fasta = merge_raw_fa(in_fa.collect(),'ont')

  // This is required to convert the ArrayList output to a Channel 
  raw_fa     = raw_fasta.fa.flatten().map { sample -> file(sample) }

  genotypes = faToGenotype("ont", raw_fa,  pubZFs, pubAlleles)

  emit:
  genotypes
  }

workflow genotypeFromMerge {
  take:
  pubZFs
  pubAlleles

  main:
  ontFA    = Channel.fromPath("${params.ontfa}/*.fa")
  pbFA     = Channel.fromPath("${params.pbfa}/*.fa")

  raw_fasta  = merge_raw_fa(ontFA.join(pbFA, remainder: true).collect(),'pb_ont')

  // This is required to convert the ArrayList output to a Channel !
  raw_fa     = raw_fasta.fa.flatten().map { sample -> file(sample) }

  genotypes = faToGenotype("final", raw_fa,  pubZFs, pubAlleles)

  emit:
  genotypes
  }

// MAIN WORKFLOW
workflow {

  aData    = get_known_human_prdm9_data()
  bsAlleles = get_jeffreys_blood_and_sperm_variants_only()
  
  gt_pacbio   = genotypeFromPacbio(aData.hsZFs,aData.hsAlleles)
  gt_nanopore = genotypeFromONT(aData.hsZFs,aData.hsAlleles)
  gt_final    = genotypeFromMerge(aData.hsZFs,aData.hsAlleles)
  
  prZFAData = mergePublishedAndFoundAlleles(gt_final,aData.hsZFs,aData.hsAlleles)

  figGels   = drawGelQuantificationFigure(gt_final)

  caveatFig = analyzeDiscordancesPB_v_ONT(gt_nanopore, gt_pacbio)

  // Infer alleles that bind common seqs
  //prZFAA    = assessZFsThatBindSimilarSequences(prZFAData.allele,prZFAData.zf)
  prZFAA    = dnaToPeptide(prZFAData.allele, prZFAData.zf, bsAlleles.alleles)
  ac        = makeACtypesPlot(prZFAA.alleles, prZFAA.zfs, bsAlleles.alleles)

  trioGTs       = checkTrios(gt_final)
  newGTs        = parseNewGTsONLY(gt_final)
  gt_validation = checkNewAllelesWithShortReads(newGTs.tab.splitText( by: 5, file: true ),
                                                prZFAA.zfs)
  figValid      = drawShortReadValidationFig(gt_validation.tab.collect())
  
  fig1Data  = drawFigure1(gt_nanopore, gt_pacbio)
  fig2Data  = drawFigure2(gt_final, prZFAA.alleles)

  // Human PrZFA network
  prZFARel  = inferRelatednessOfAlleles(prZFAData.names.splitText( by: 8, file: true ),prZFAA.alleles)
  
  // merge relatedness files
  allRelatednessData = prZFARel.pub.collectFile(name: 'PrZFA_relatedness_ALLhuman.tab', newLine: true)
  
  fig3Data  = drawRelatednessFigure(prZFARel.pub.collect(), 
                                    prZFAData.allele,
                                    prZFAA.alleles,
                                    bsAlleles.alleles,
                                    gt_final)

  // Analysis of linked haplotypes
  aLink     = getLinkedAlleles(gt_final,prZFAA.alleles)
  figAssoc  = drawAssociationsFigure(aLink.prdm9.collect(),
                                    aLink.cmh.collect(),
                                    aLink.type.collect(),
                                    aLink.clust.collect(),
                                    aLink.map.collect(),
                                    aLink.snpseq.collect())

  // Hotspots figure : SSDS & Motifs
  hotspots   = Channel.fromPath("${params.hs}/*.bedgraph")

  hsData     = processHumanHotspots(hotspots.collect())

  memeMotifs = getMotifsMEME(hsData.memeFA.flatten(), hsData.wideFA.collect())

  fig5Data   = drawHotspotsFigure(memeMotifs.meme.collect(), hsData.table)
}