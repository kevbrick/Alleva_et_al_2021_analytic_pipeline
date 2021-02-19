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
params.outdir         = ""
params.pbfa           = ""
params.ontfa          = ""
params.accessorydir   = ""

params.genomefa       = ""
params.genomeidx      = ""
params.npeaks         = 1000
params.bigwinSize     = 1600
params.centerSz       = 300

genomeFA  = file("${params.genomefa}")
genomeIDX = file("${params.genomefa}.fai")

//log.info
log.info " "
log.info "=========================================================================="
log.info "PRDM9 genotyping PIPELINE (Version 1.0.1)                                "
log.info "=========================================================================="
log.info " "
log.info "USAGE: "
log.info " "
log.info "------------------------------------------------------------------------- "
log.info "nextflow run genotype_prdm9.nf "
log.info " --pbfa          ${params.pbfa} "
log.info " --ontfa (ONT)   ${params.ontfa} "
log.info " --outdir        ${params.outdir} "
log.info " --accessorydir  ${params.accessorydir} "
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
  //
  // ## Split FASTA files to equal sized subfolders
  // i=1
  // while read l; do
  //   mkdir -p faSplit_\$i
  //   ln -s \$l faSplit_\$((i++))
  // done< <(find -L `pwd` -name '*.pb_ont.raw.fa' | xargs -n 2)

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

  echo -e `grep allele_1_code allhap.txt`"\\tnseqs\\tnzfseqs" |perl -pi -e 's/\\s+(\\S)/\\t\$1/g' >prdm9_haplotypes.${type}.tab

  for f in *.fa; do

    if [[ ! \$f =~ "ZFs.fa" ]]; then
      nm=`basename \$f |perl -pi -e 's/^(\\S+?)_.+\$/\$1/' 2>/dev/null`

      if [[ `find -name "\$nm*ZFs.fa"` ]]; then
        sz=`grep -P "^\\>" \$f |wc -l`
        zfsz=`cat \$nm*ZFs.fa |grep \\> |wc -l`

        hap=`grep \$nm allhap.txt`

        echo -e "\$hap\\t\$sz\\t\$zfsz" >>prdm9_haplotypes.${type}.tab
        echo "Processed \$nm ..."
      else
        echo "Skipping \$nm: No ZFs found ..."
      fi
    fi
  done

  touch newAlleles.txt
  touch newZFs.txt
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

process getMousePrZFAs{
  publishDir "${params.outdir}/annotation/mousePrZFA", mode: 'copy', overwrite: true, pattern: 'PrZFA*txt'

  output:
  path('PrZFA_allele_codes.mouse.txt', emit: allele)
  path('PrZFA_ZF_sequences.mouse.txt', emit: zf)
  path('distinct*txt',                 emit: names)

  """
  ln -s ${params.accessorydir} accessoryFiles

  cp accessoryFiles/otherdata/mouse_PrZFA_alleles_in_Pubs.txt .

  perl accessoryFiles/scripts/getMMPrZFASequences.pl

  perl accessoryFiles/scripts/convertFAToMousePrZFAs.pl

  cut -f1 PrZFA_allele_codes.mouse.txt |uniq >distinct_PrZFA_allele_names_mouse.txt

  """
  }

process drawFigure1{
  publishDir "${params.outdir}/figures", mode: 'copy', overwrite: true

  input:
  path(bonito_dets)
  path(pacbio_dets)

  output:
  path('*png', emit: png)
  path('*pdf', emit: pdf)

  """
  ln -s ${params.accessorydir} accessoryFiles

  R --no-save <${params.accessorydir}/scripts/drawAllevaFigure1.R
  rm -f Rplots.pdf
  """
  }

process drawFigure2{
  publishDir "${params.outdir}/figures", mode: 'copy', overwrite: true

  input:
  path(dets)

  output:
  path('Alleva*png', emit: png)
  path('Alleva*pdf', emit: pdf)

  """
  ln -s ${params.accessorydir} accessoryFiles
  cp accessoryFiles/scripts/genericFunctions.R .

  R --no-save <${params.accessorydir}/scripts/drawAllevaFigure2.R
  """
  }

process drawFigure3{
  publishDir "${params.outdir}/figures",    mode: 'copy', overwrite: true, pattern: '*png'
  publishDir "${params.outdir}/figures",    mode: 'copy', overwrite: true, pattern: '*pdf'
  publishDir "${params.outdir}/annotation", mode: 'copy', overwrite: true, pattern: '*txt'

  input:
  path(rtable)
  path(alleles)

  output:
  path('Alleva_et_al_Fig*p*',         emit: mainFig)
  path('Alleva_et_al_Sup*p*',         emit: suppFigs)
  path('PrZFA_relatedness.human.tab', emit: rtab)

  """
  ln -s ${params.accessorydir} accessoryFiles

  grep event PrZFA_relatedness_ALLhuman.tab       |head -n1                                        >PrZFA_relatedness.human.tab
  sort -k1,1 -k2,2 PrZFA_relatedness_ALLhuman.tab |grep -v event |perl -lane 'print \$_ if (\$_)' >>PrZFA_relatedness.human.tab

  grep -P 'TRUE\\s*\$' PrZFA_alleles.txt >PrZFA_alleles.foundInPops.txt

  cp accessoryFiles/scripts/genericFunctions.R . 

  R --no-save <${params.accessorydir}/scripts/drawAllevaFigure3.R
  rm -f Rplots.pdf
  """
  }

process drawFigure3MMSupp{
  publishDir "${params.outdir}/figures",    mode: 'copy', overwrite: true, pattern: '*png'
  publishDir "${params.outdir}/figures",    mode: 'copy', overwrite: true, pattern: '*pdf'
  publishDir "${params.outdir}/annotation", mode: 'copy', overwrite: true, pattern: '*txt'

  input:
  path(rtable)
  path(alleles)

  output:
  path('Allev*.p??',                  emit: mmFig)
  path('PrZFA_relatedness.mouse.tab', emit: rtab)

  """
  ln -s ${params.accessorydir} accessoryFiles

  grep event PrZFA_relatedness_ALLmouse.tab       |head -n1                                        >PrZFA_relatedness.mouse.tab
  sort -k1,1 -k2,2 PrZFA_relatedness_ALLmouse.tab |grep -v event |perl -lane 'print \$_ if (\$_)' >>PrZFA_relatedness.mouse.tab

  cp ${alleles} PrZFA_alleles.foundInPops.txt

  perl -lane '\$F[0] =~ /^(\\S+?)_(\\S+)\$/;
             (\$str,\$id) = (\$1,\$2);
             if (\$id =~ s/::(\\S+)\$//){\$pub = \$1}else{\$pub = "NA"};
             print join("\\t",\$str,\$F[0],\$id,\$pub,\$F[1])' PrZFA_allele_codes.mouse.txt >PrZFA_allele_details.mouse.txt

  R --no-save <${params.accessorydir}/scripts/drawAllevaFigure3_MouseNetworkSupp.R
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

  perl -lane 'print join("\\t",\$F[4],\$F[6]);
              print join("\\t",\$F[5],\$F[7])' ${dets} |sort -V |grep -v Unk |grep -v NA |grep -v allele |uniq >PrZFA_alleles.foundInPops.txt

  sort -V ${pubZFs}     PrZFA_ZFs.foundInPops.txt     |uniq >PrZFA_ZFs.txt
  sort -V ${pubAlleles} PrZFA_alleles.foundInPops.txt |uniq >PrZFA_alleles.tmp

  ## Add logical fields for published / pop / new alleles
  grep -f    ${pubAlleles} PrZFA_alleles.tmp |perl -lane 'chomp; print "\$_\\tTRUE"'  >aPub.txt
  grep -v -f ${pubAlleles} PrZFA_alleles.tmp |perl -lane 'chomp; print "\$_\\tFALSE"' >aNoPub.txt
  sort -V aPub.txt aNoPub.txt >PrZFA_alleles.vPub.tmp

  grep -f    PrZFA_alleles.foundInPops.txt PrZFA_alleles.vPub.tmp |perl -lane 'chomp; print "\$_\\tTRUE\\t"'  >aPop.txt
  grep -v -f PrZFA_alleles.foundInPops.txt PrZFA_alleles.vPub.tmp |perl -lane 'chomp; print "\$_\\tFALSE\\t"' >aNoPop.txt
  sort -V aPop.txt aNoPop.txt >PrZFA_alleles.txt

  cut -f1 PrZFA_alleles.txt |uniq >distinct_PrZFA_allele_names.txt

  """
  }

process inferRelatednessOfAlleles{
  input:
  path(alleles_to_check)
  path(PrZFA_alleles)
  path(PrZFA_ZFs)

  output:
  path('PrZFA_relatedness*tab', emit: pub)

  """
  ln -s ${params.accessorydir} accessoryFiles

  perl ${params.accessorydir}/scripts/checkPrdm9Recombinants.pl --a ${PrZFA_alleles} --z ${PrZFA_ZFs} --c ${alleles_to_check} --o PrZFA_relatedness.\$RANDOM\$RANDOM.tab
  """
  }

process inferRelatednessOfMouseAlleles{
  input:
  path(alleles_to_check)
  path(PrZFA_alleles)
  path(PrZFA_ZFs)

  output:
  path('PrZFA_relatedness*tab', emit: pub)

  """
  ln -s ${params.accessorydir} accessoryFiles

  perl ${params.accessorydir}/scripts/checkPrdm9Recombinants.pl --a ${PrZFA_alleles} --z ${PrZFA_ZFs} --c ${alleles_to_check} --o PrZFA_relatedness.\$RANDOM\$RANDOM.tab
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

process getMotifsChipmunk {

  publishDir "${params.outdir}/humanHS/motifs/chipmunk",   mode: 'copy', overwrite: true, pattern: '*.chipmunk*'
  publishDir "${params.outdir}/humanHS/motifs/chipmunk",   mode: 'copy', overwrite: true, pattern: '*.dichipmunk*'
  //publishDir "${params.outdir}/humanHS/motifs/chiphorde",   mode: 'copy', overwrite: true, pattern: '*.chiphorde*'
  //publishDir "${params.outdir}/humanHS/motifs/dichiphorde", mode: 'copy', overwrite: true, pattern: '*.dichiphorde*'
  publishDir "${params.outdir}/humanHS/motifs",            mode: 'copy', overwrite: true, pattern: '*png'
  publishDir "${params.outdir}/humanHS/motifs",            mode: 'copy', overwrite: true, pattern: '*density.txt'

  tag { fasta }

  input:
  path(fasta)
  path(wideFA)

  output:
  path("*.chipmunk*",    optional: true, emit: cm)
  path("*.dichipmunk*",  optional: true, emit: dcm)
  //path("*.chiphorde*",   optional: true, emit: ch)
  //path("*.dichiphorde*", optional: true, emit: dch)
  path("*motifdensity*", optional: true, emit: motifdist)

  script:
  def name  = fasta.name.replaceFirst(".chipmunk.fa","")

  """
  ruby /opt/chipmunk/run_chipmunk8.rb     ${name}.chipmunk    24 12 yes 1.0 m:${fasta} || true
  ruby /opt/chipmunk/run_dichipmunk8.rb   ${name}.dichipmunk  24 12 yes 1.0 m:${fasta} || true
  #ruby /opt/chipmunk/run_chiphorde8.rb    ${name}.chiphorde   24:12,24:12 filter verbose=yes 1.0 m:${fasta} || true
  #ruby /opt/chipmunk/run_dichiphorde8.rb  ${name}.dichiphorde 24:12,24:12 filter verbose=yes 1.0 m:${fasta} || true

  for pcm in *.pcm; do
    png=\${pcm/pcm/FWD.png}
    ruby /opt/chipmunk/pmflogo3.rb \$pcm \$png 1 100 400 weblogo no
    png=\${pcm/pcm/REV.png}
    ruby /opt/chipmunk/pmflogo3.rb \$pcm \$png 1 100 400 weblogo yes
  done

  for pcm in *.dpcm; do
    png=\${pcm/dpcm/FWD.png}
    ruby /opt/chipmunk/dpmflogo3.rb \$pcm \$png 100 400
    png=\${pcm/dpcm/REV.png}
    ruby /opt/chipmunk/dpmflogo3_rc.rb \$pcm \$png 100 400
  done

  mfile=${name}".chipmunkmotifdensity.txt";
  echo -e "hotspots\\ttype\\tmotifsrc\\tscore\\tpos\\tstrand\\tcs\\tfrom\\tto\\tstrength\\tallele" >\$mfile

  for fa in *.bigwin.fa; do
    genotype=\${fa/.bigwin.fa};

    for pwm in *.pwm; do
      motifsource=\${pwm/.chipmunk.pwm/}
      java -cp \$SARUSJAR ru.autosome.SARUS \$fa \$pwm 1 transpose |perl -lane 'if (\$_ =~ /\\>\\s*(.+)\$/){@nm = split(":",\$1); next}; print join("\\t","'\$genotype'","mono","'\$motifsource'",\$_,@nm)' >>\$mfile
    done

    for dpwm in *.dpwm; do
      motifsource=\${pwm/.chipmunk.dpwm/}
      java -cp \$SARUSJAR ru.autosome.di.SARUS \$fa \$dpwm 1 transpose |perl -lane 'if (\$_ =~ /\\>\\s*(.+)\$/){@nm = split(":",\$1); next}; print join("\\t","'\$genotype'","di","'\$motifsource'",\$_,@nm)' >>\$mfile
    done
  done

  """
  }

process drawFigure4 {

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
  rm -f Rplots.pdf
  """
  }

//#############################################
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

workflow genotypeFromBonito {
  take:
  pubZFs
  pubAlleles

  main:
  bonito_fa = Channel.fromPath("${params.ontfa}/*.fa")
  raw_fasta = merge_raw_fa(bonito_fa.collect(),'ont')

  // This is required to convert the ArrayList output to a Channel !
  raw_fa     = raw_fasta.fa.flatten().map { sample -> file(sample) }

  genotypes = faToGenotype("bonito", raw_fa,  pubZFs, pubAlleles)
  //genotypes = faToGenotype("bonito", raw_fasta,  pubZFs, pubAlleles)

  emit:
  genotypes
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

workflow genotypeFromGuppy {
  take:
  pubZFs
  pubAlleles

  main:
  guppy_fa  = Channel.fromPath("${params.guppyfa}/1D2/*.fa")
  raw_fasta = merge_raw_fa(guppy_fa.collect())

  // This is required to convert the ArrayList output to a Channel !
  raw_fa     = raw_fasta.fa.flatten().map { sample -> file(sample) }

  genotypes = faToGenotype("guppy", raw_fas,  pubZFs, pubAlleles)

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

  gt_bonito = genotypeFromBonito(aData.hsZFs,aData.hsAlleles)
  gt_pacbio = genotypeFromPacbio(aData.hsZFs,aData.hsAlleles)
  //gt_guppy  = genotypeFromGuppy(aData.hsZFs,aData.hsAlleles)
  gt_final  = genotypeFromMerge(aData.hsZFs,aData.hsAlleles)

  prZFAData = mergePublishedAndFoundAlleles(gt_final,aData.hsZFs,aData.hsAlleles)
  mmZFAData = getMousePrZFAs()

  figGels   = drawGelQuantificationFigure(gt_final)

  fig1Data  = drawFigure1(gt_bonito, gt_pacbio)
  fig2Data  = drawFigure2(gt_final)

  // Human PrZFA network
  prZFARel  = inferRelatednessOfAlleles(prZFAData.names.splitText( by: 4, file: true ),prZFAData.allele,prZFAData.zf)
  fig3Data  = drawFigure3(prZFARel.pub.collectFile(name: 'PrZFA_relatedness_ALLhuman.tab', newLine: true), prZFAData.allele)

  // Mouse PrZFA network
  mmZFARel  = inferRelatednessOfMouseAlleles(mmZFAData.names.splitText( by: 4, file: true ),mmZFAData.allele,mmZFAData.zf)
  fig3MMSu  = drawFigure3MMSupp(mmZFARel.pub.collectFile(name: 'PrZFA_relatedness_ALLmouse.tab', newLine: true), mmZFAData.allele)

  // Figure 4 : SSDS & Motifs
  hotspots   = Channel.fromPath("${params.hs}/*.bedgraph")
  //hapTable   = Channel.fromPath("${params.hapdir}/haplotypesTable/prdm9_haplotypes.tab")
  pubData    = Channel.fromPath("${params.hapdir}/publishedData/*")

  hsData     = processHumanHotspots(hotspots.collect())

  memeMotifs = getMotifsMEME(hsData.memeFA.flatten(), hsData.wideFA.collect())

  fig4Data   = drawFigure4(memeMotifs.meme.collect(), hsData.table)
  //cmMotifs   = getMotifsChipmunk(hsData.chipmunkFA.flatten(), hsData.wideFA.collect())
}


// process getErrorRates{
//   publishDir "${params.outdir}/figures", mode: 'copy', overwrite: true
//
//   input:
//   path(guppy)
//   path(bonito)
//   path(pacbio)
//
//   output:
//   path('*png', emit: png, optional: true)
//   path('*pdf', emit: pdf, optional: true)
//
//   """
//   cut  -f1       prdm9_haplotypes.guppy.tab           >g1.txt
//   grep -f g1.txt prdm9_haplotypes.bonito.tab |cut -f1 >g2.txt
//   grep -f g2.txt prdm9_haplotypes.pacbio.tab |cut -f1 >g3.txt
//   grep -f g3.txt prdm9_haplotypes.final.tab  |grep -P 'A\/A' |cut -f1 |shuf |head -n20 >r20.txt
//
//   """
//   }
