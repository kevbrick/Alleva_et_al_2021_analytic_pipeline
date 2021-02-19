# Get an EMBOSS factory
use strict;
use Bio::Seq;
use Bio::SeqIO;
use Bio::AlignIO;
use Bio::SimpleAlign;
use Bio::LocatableSeq;
use Math::Round qw (round);
use List::UtilsBy qw(max_by);

use Getopt::Long;


GetOptions ('fa=s'     => \(my $seq_input_fa),
	    'id=s'     => \(my $sample_id),
	    'c=s'      => \(my $consensus_threshold = 0.7),
	    'n=i'      => \(my $n_seqs_required = 20));

die("FASTA ($seq_input_fa) does not exist!") unless (-e $seq_input_fa);
die("Invalid consensus threshold (0 - 1)")   if ($consensus_threshold < 0 || $consensus_threshold >1);
die("Sample name required (--id)!")          unless ($sample_id);
die("Cannot find R (load module)!")          unless (`which R`);
die("Cannot find replaceUnicode script in current folder!")          unless (-e 'replaceUnicode');
#die("Cannot find rhs.fa in current folder!")          unless (-e 'rhs.fa');
#die("Cannot find lhs.fa in current folder!")          unless (-e 'lhs.fa');
#die("Cannot find zfA.fa in current folder!")          unless (-e 'zfA.fa');
#die("Cannot find zf.fa in current folder!")           unless (-e 'zf.fa');

my ($gap_open_penalty,$gap_extend_penalty) = (20,2);
my $rand_name_base                         = "zfMatch_".int(rand()*1000000000000000);

#my ($seq_input_fa,$sample_id,$percent_required,) = @ARGV;

## make BioPerl object to load ZF array $fasta
## This is simply a fasta made from column 10 of the ccs BAM
my $ccs_fasta                              = Bio::SeqIO->new(-format => 'fasta',
                                                             -file   => $seq_input_fa);

open OUTZFS, '>', "prdm9_alleles_from_ccs.".$sample_id.".raw";
open OUTHAP, '>', "prdm9_alleles_from_ccs.".$sample_id.".haplotypes.txt";
open LOG, '>', "prdm9_alleles_from_ccs.log";

my ($totSeqs, $noZFSeqs, $goodSeqs);
my ($fasta_seq_count, %zf_hash, $ok_sequence_count, %zf_array_count_by_size);
my ($zf_array_count_total, %consensus_zfs, %two_alleles, %zf_array_sizes_to_use, %zf_aln);
my (@genotypes,$got_diploid_genotype,$genotyping_problem);

## This loop will find the ZFs in each sequence, check that the flanks are found, and build a
## diploid consensus sequence for each. We loop over each sequence in the fasta
while ( my $one_ccs_seq = $ccs_fasta->next_seq ) {

  ## Unnecessary precaution, but I'll leave it here.
  last if (++$fasta_seq_count > 500000);

  $totSeqs++;

  ## Triage reads with no ZFs
  if ($one_ccs_seq->seq !~ /(TGTGGGC|GCCCACA|GGGCTTT|AAAGCCC|GGGGAGA|TCTCCCC)/){
    $noZFSeqs++;
    next;
  }

  my $seq2use = $one_ccs_seq->seq;
  if ($one_ccs_seq->seq =~ /(GCCCACA|AAAGCCC|TCTCCCC)/){
    $seq2use = $one_ccs_seq->seq ;
  }

  ## Calls the get_zf_array function. This will find sequences that have a
  ## complete zf array and return details of the array
  my ($number_of_zfs, $both_edges_found, $seq_of_entire_zf_array, @zfs) = getZFArray($seq2use, $one_ccs_seq->id);

  next unless ($number_of_zfs);

  my ($zf_array_sizes_known);

  ## Select for COMPLETE ZF arrays - ONLY proceed if the ccs contains BOTH the left and right flanking sequences.
  if ($both_edges_found){

    ## We're going to loop through the ZFs, so set counter to zero
    my $current_zf_pos = 0;

    ## This keeps track of the number of zf arrays by size
    $zf_array_count_by_size{$number_of_zfs}++;

    $zf_array_count_total++;
    $goodSeqs++;

    ## Try to determine the likely sizes of the two alleles.
    $zf_array_sizes_known = check_zf_array_sizes($zf_array_count_total, \%zf_array_count_by_size, \%two_alleles, \%zf_array_sizes_to_use);

    ## If we know the zf sizes, then stop processing zf arrays that are the wrong size
    next unless ($zf_array_sizes_to_use{$number_of_zfs} || not ($zf_array_sizes_known));

    ## Next, what we do is keep track of the ZF sequences found at each position in the array
    ## We will continue until we reach a reasonable $consensus

    ## Create bioperl alignment if necessary
    $zf_aln{$number_of_zfs} = Bio::SimpleAlign->new() unless ($zf_aln{$number_of_zfs});

    my $seq = Bio::LocatableSeq->new(-seq   => $seq_of_entire_zf_array,
                                     -id    => "zfa:".$number_of_zfs.":".$zf_array_count_total,
                                     -start => 1,
                                     -end   => length($seq_of_entire_zf_array));

    $zf_aln{$number_of_zfs}->add_seq($seq);
  }
  if ($goodSeqs >= $n_seqs_required && !($goodSeqs % 5)){
    if ($zf_array_sizes_known){
      ($genotyping_problem,@genotypes) = try_to_guess_genotypes(\%zf_aln,\%zf_array_sizes_to_use);
      if ($genotyping_problem){
        print STDERR "********************************************************\n";
        print STDERR "* Loop : $zf_array_count_total\n";
        print STDERR "* Genotyping problem: $genotyping_problem\n";
        print STDERR "********************************************************\n";
      }else{
        for my $gt(@genotypes){
          my ($zf_array_sequence,$totZF);
          for my $nzf(sort {$a <=> $b} keys(%{$gt})) {
            print OUTZFS join("\t",$nzf,$gt->{$nzf})."\n";
            $zf_array_sequence .= $gt->{$nzf};
            $totZF++;
          }
          print OUTHAP join("\t",$sample_id,$totZF,($#genotypes+1),$zf_array_sequence)."\n";
        }
        close OUTZFS;
        close OUTHAP;
        exit;
      }
    }else{
      print STDERR "ZF Array size ambiguity ... continuing ... \n";
    }
  }
}

close LOG;

sub getZFArray{
  my ($inSeq,$rname) = @_;

  my $tmpBase  = "/lscratch/".$ENV{'SLURM_JOBID'}."/$rand_name_base";
  my $outWater = "$tmpBase.water";
  my $inFasta  = "$tmpBase.input.fa";

  open FA, '>', $inFasta;
  print FA ">prdm9\n$inSeq";
  close FA;

  makeFASTAs("$tmpBase");

  my ($sSeq,$aScore,$aFrom,$aTo) = getPosOnRead("$tmpBase.ZFA.fa",$inFasta);
  my ($lSeq,$lScore,$lFrom,$lTo) = getPosOnRead("$tmpBase.LHS.fa",$inFasta);
  my ($rSeq,$rScore,$rFrom,$rTo) = getPosOnRead("$tmpBase.RHS.fa",$inFasta);

  ##print STDERR "Aligning Read $fasta_seq_count ...\n";
  system("matcher -alt 30  -gapopen $gap_open_penalty -gapextend $gap_extend_penalty -asequence $tmpBase.ZF.fa -bsequence $inFasta -outfile $outWater -aformat pair 2>/dev/null");

  my @okArray;
  # Now you might want to get the alignment
  my $in = Bio::AlignIO->new(-format => 'emboss',
                             -file   => $outWater);

  my ($zfSeq, %zfs, $outSeq,$skip,$padleft,$padright,$nAln);
  my ($arrLOK, $arrROK, $both_edges_found, $zfarray_seq, $num_zfs);

  my ($minPos,$maxPos) = (1e9,-1e9);
  while ( my $aln = $in->next_aln ) {
    ($skip,$padleft,$padright,$nAln) = (0,"","",0);
    #if ($aln->length == 84){
    my (@qry,@tgt,$s,$f,$t);
    ##print STDERR "Parsing Alingment # ".(++$nAln)." ; Read $fasta_seq_count ...\n";
    foreach my $seq ($aln->each_seq) {
      if ($seq->id eq 'prdm9'){

        $f = $seq->start;
        $t = $seq->end;
        $s =  $seq->seq;
        #$s =  $seq->subseq($f, $t);
        #print STDERR join("\t",$s,$f,$t,$aln->score."\n");

        $zfs{$f} -> {from}  = $f;
        $zfs{$f} -> {to}    = $t;
        $zfs{$f}->{skip}    = 0;
        $zfs{$f} -> {score} = $aln->score;

        #$zfSeq .= $s;
        $minPos = ($minPos < $f)?$minPos:$f;
        $maxPos = ($maxPos > $f)?$maxPos:$f;

        if ($zfs{$f} -> {from} < ($lTo-5)){
          $zfs{$f}->{skip}    = 1 ;
          my $rere = 1;
        }

        #for my $i(($f-$aFrom)..($t-$aFrom)){$okArray[$i] = 1};
        @tgt = split("",$s);
      }

      if ($seq->id eq 'zf'){
        #$skip++ unless ($seq->start == 1 && $seq->end == 84);
        @qry = split("",$seq->seq);
        #$skip++ if ($seq->end != 84 || $seq->start != 1);
        if ($seq->start != 1){
          $padleft  = "N" x ($seq->start-1) ;
        }
        if ($seq->end != 84){
          $padright = "N" x (84-$seq->end) ;
          my $rrrr = 1;
        }
      }

      #$zfs{$f}->{skip} = 1 if ($f && $skip);

    }

    my $mySeq = $padleft;

    for my $i(0..$#qry){
      $mySeq .= $tgt[$i] if ($qry[$i] ne "-");
    }

    #return 0 if ($skip);
    $zfs{$f}->{seq}  = $mySeq.$padright;

    if (length($zfs{$f}->{seq}) < 83 || length($zfs{$f}->{seq}) >84 ){
      my $rasd = 1;
    }
    #print STDERR join("",@qry)."\n";
    #print STDERR join("",@tgt)."\n";
    #print STDERR $zfs{$f}->{from}." : ".$mySeq." : ".$zfs{$f}->{to}."\n";
    #print STDERR ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";

  }

  my $kPrev = 0;
  my @outDetail;

  my (%zfOK,%numnum,@ZFnumz);
  foreach my $k (sort {$a <=> $b} (keys %zfs)) {
    my $nZF = round(($zfs{$k}->{from}-($lTo+1))/84)+1;

    my $skipZF;
    $skipZF++ if ($zfOK{$nZF} && $zfOK{$nZF}->{score} > $zfs{$k}->{score});

    unless ($skipZF){#} || $zfs{$k}->{score} < 200){
      $zfOK{$nZF}->{seq}     = $zfs{$k}->{seq};
      $zfOK{$nZF}->{zf}      = round(($zfs{$k}->{from}-$minPos)/84);
      $zfOK{$nZF}->{zfaFrom} = $zfs{$k}->{from}-$minPos;
      $zfOK{$nZF}->{zfaTo}   = $zfs{$k}->{to}-$minPos;
      $zfOK{$nZF}->{from}    = $zfs{$k}->{from};
      $zfOK{$nZF}->{to}      = $zfs{$k}->{to};
      $zfOK{$nZF}->{score}   = $zfs{$k}->{score};
      $zfOK{$nZF}->{skip}    = $zfs{$k}->{skip};
      $zfOK{$nZF}->{lOK}     = ($lTo   == ($zfs{$k}->{from}-1))?1:0;
      $zfOK{$nZF}->{rOK}     = ($rFrom == ($zfs{$k}->{to}+1))?1:0;
      push @ZFnumz, $nZF unless ($numnum{$nZF}++);
      $arrLOK++ if ($zfOK{$nZF}->{lOK});
      $arrROK++ if ($zfOK{$nZF}->{rOK});
      $both_edges_found = ($arrLOK && $arrROK)?1:0;
      $zfarray_seq .= $zfs{$k}->{seq};
    }

  }

  return unless($both_edges_found);

  $zfarray_seq =~ s/\-/X/g;

  #print STDERR "Read # $fasta_seq_count :\n";
  #foreach my $k (sort {$a <=> $b} (keys %zfs)) {
  my ($lChk,$rChk,$allZFchk,@zfSequences);
  foreach my $k (@ZFnumz) {
    #if ($zfOK{$k}->{skip} eq 'OK'){
        #return 0 if ($kPrev && (($zfs{$k}->{from} - 1) == $zfs{$k}->{to}));
        $kPrev = $k;

        next if ($rChk);

        $lChk++ if ($zfOK{$k}->{lOK});
        $rChk++ if ($zfOK{$k}->{rOK});

        if ($lChk){
          if ($allZFchk){
            $allZFchk = ($allZFchk == $k-1)?$k:999;
          }else{
            $allZFchk = $k;
          }
        }

        if ($lChk){
          push @outDetail, join("\t",$fasta_seq_count,
                                    $zfOK{$k}->{seq},
                                    $k,
                                    $zfOK{$k}->{zfaFrom},
                                    $zfOK{$k}->{zfaTo},
                                    $lTo." -- ".$rFrom,
                                    $zfOK{$k}->{from},
                                    $zfOK{$k}->{to},
                                    $zfOK{$k}->{score},
                                    $zfOK{$k}->{skip}?"SKIP":($zfOK{$k}->{score}<200?"LOSCORE":"OK"),
                                    $zfOK{$k}->{lOK}?"LOK":"noL",
                                    $zfOK{$k}->{rOK}?"ROK":"noR");
          push @zfSequences, $zfOK{$k}->{seq};

          $num_zfs++;
          $zfSeq .= $zfOK{$k}->{seq};
        }
      #}
  }

  for my $o(@outDetail){print STDERR $o."\t".join("\t",$totSeqs,$noZFSeqs,$goodSeqs,"SO")."\n"} ;# if ($allZFchk < 999 && $rChk)}
  #print STDERR "---------------------------------------------------------------------------------------------------\n";
  #print ">".$rname."\n".to84($zfSeq);

  $zfarray_seq = join("",@zfSequences);
  $zfarray_seq =~ s/\-/X/g;

  return($num_zfs, $both_edges_found, $zfarray_seq, @zfSequences);
}

sub to84{
  my $inLine = shift;
  my $lnRet;
  while ($inLine =~ s/^(.{84})//){
    $lnRet .= $1."\n";
  }
  $lnRet .= $1."\n";
  return($lnRet);
}

sub getPosOnRead{
  my ($checkFA,$iFA) = @_;

  my $tf = "/lscratch/".$ENV{'SLURM_JOBID'}."/$rand_name_base.out";
  ## ALIGN TO GET "A" ALLELE FIRST
  system("matcher -alt 30 -gapopen $gap_open_penalty -gapextend $gap_extend_penalty -asequence $checkFA -bsequence $iFA -outfile $tf -aformat pair 2>/dev/null");

  # Now you might want to get the alignment
  my $in = Bio::AlignIO->new(-format => 'emboss',
                             -file   => $tf);

  my ($iSeq,$iScore,$iFrom,$iTo);

  while ( my $aln = $in->next_aln ) {
      foreach my $seq ($aln->each_seq) {
        if ($seq->id eq 'prdm9'){
          $iFrom = $seq->start;
          $iTo = $seq->end;
        }
        return("",999,$iFrom,$iTo) if ($iFrom);
      }
  }
  return(0,0,0,0);
}

sub check_zf_array_sizes{
  my ($tot,$sz,$two_X_alleles,$ok_sizes,$zf_known) = @_;

  return 0 if ($tot < $n_seqs_required);

  my @zf_array_rank = sort { $$sz{$b} <=> $$sz{$a} } keys(%{$sz});

  if ($$sz{$zf_array_rank[0]} > $tot*$consensus_threshold){
    $$ok_sizes{$zf_array_rank[0]}++;
    $$two_X_alleles{$zf_array_rank[0]}++;
    return 1;
  }

  ## If the top two ZFs represent N% of all tested ccs
  if (($$sz{$zf_array_rank[0]} + $$sz{$zf_array_rank[1]}) > $tot*$consensus_threshold){
    $$ok_sizes{$zf_array_rank[0]}++;
    $$ok_sizes{$zf_array_rank[1]}++;
    return 1;
  }

  ## If we got here, we still don't know
  return 0;
}

sub do_we_have_a_consensus_diploid_sequence{
  my ($h_consensus, $h_sizes) = @_;

  return 0 unless (%{$h_consensus});

  my ($num_ok, $num_required);

  for my $sz(keys(%{$h_sizes})){
    my @zf_pos_with_consensus = keys(%{$$h_consensus{$sz}});
    $num_ok++ if (($#zf_pos_with_consensus+1) == $sz);
    $num_required++;
  }

  return ($num_ok && ($num_ok == $num_required))?1:0;
}

sub try_to_guess_genotypes{

  my ($aln_all, $array_sizes) = @_;

  my ($problems, @returnGenotypes);

  my @n_asizes = keys(%{$array_sizes});
  my $minimum_alleles_required = ($#n_asizes+1);

  for my $array_size(keys(%{$array_sizes})){
    my $aln = $$aln_all{$array_size};

    my (@seq_array, @locSeqArray);

    for my $s($aln->each_seq){
      push @seq_array, Bio::Seq->new(-seq => $s->seq,-id => $s->id);
      push @locSeqArray, $s;
    }

    my @heterozygosity_pos;
    my $pos = 0;

    for my $pc_id($aln->consensus_conservation()){
      if ($pc_id < 80){
        #print $pos."\n" ;
        push @heterozygosity_pos, $pos;
      }
      $pos++;
    }

    ## No reason to think its a het ... return consensus.
    if (@heterozygosity_pos){
      my $short_aln = Bio::SimpleAlign->new();
      my (@seqArrayShort, @ss);
      my $s = 0;

      for my $sq(@seq_array){
        my $selectSeq;
        $selectSeq .= $sq->subseq($_+1,$_+1) foreach(@heterozygosity_pos);
        push @ss, $selectSeq;

        my $seq = Bio::LocatableSeq->new(-seq => "$selectSeq",
                       -id  => "selseq".$s++,
                       -alphabet => 'dna',
                       -start => 1,
                       -end   => length($selectSeq));

        $short_aln->add_seq($seq);

        my $sSeq = Bio::Seq->new(-seq => "$selectSeq",
                                 -id  => "seq".$s++);
        push @seqArrayShort, $sSeq;
      }

      my $matrixFile = "/lscratch/".$ENV{SLURM_JOBID}.'/tst_'.int(rand()*1000000000000).'.csv';
      clusterMe($matrixFile,@ss);
      system("Rscript getClustersFromDistMatrix.R --distmatrix $matrixFile 2>/dev/null");

      my $clusterFile = $matrixFile.'.clusters.csv';
      open IN, $clusterFile;
      open my $CF, '<', $clusterFile;
      chomp(my @clusters = <$CF>);
      close $CF;

      my %clustAln;
      $clustAln{1} = Bio::SimpleAlign->new();
      $clustAln{2} = Bio::SimpleAlign->new();

      my $nSeq;
      for my $i(@clusters){
          $clustAln{$i}->add_seq($locSeqArray[$nSeq++])
      }

      my (%consensus, %cons1, %cons2, $ambiguity1, $ambiguity2);
      $ambiguity1 = splitToZFs($clustAln{1}, \%cons1);
      $ambiguity2 = splitToZFs($clustAln{2}, \%cons2);

      $problems .= ',ambiguity_in_haplotype1'  if ($ambiguity1);
      $problems .= ',ambiguity_in_haplotype2'  if ($ambiguity2);

      push @returnGenotypes, \%cons1;
      push @returnGenotypes, \%cons2;
    }else{
      my %consensus;
      splitToZFs($aln, \%consensus);

      push @returnGenotypes, \%consensus;
    }
  }

  $problems .= ',too_many_genotypes('.($#returnGenotypes+1).')' if ($#returnGenotypes > 1); ## Zero-based
  $problems .= ',too_few_genotypes('.($#returnGenotypes+1).')'  if (($#returnGenotypes+1) < $minimum_alleles_required );

  return($problems,@returnGenotypes);
}

###############################################
sub clusterMe{
  my ($matrix_csv,@s) = @_;
  open OUT, '>', $matrix_csv;
  my $matrix;
  for my $ii(0..$#s){
    for my $jj(0..$#s){
      my $dist = 0;
      if ($s[$ii] eq $s[$jj]){
        $matrix .= "0,";
      }else{
        my @seq1 = split(//,$s[$ii]);
        my @seq2 = split(//,$s[$jj]);
        for my $kk(0..$#seq1){
          $dist++ if ($seq1[$kk] ne $seq2[$kk]);
        }

        $matrix .= "$dist,";
      }
    }
    $matrix =~ s/[,\s]+$//;
    print OUT $matrix."\n";
    $matrix="";
  }

  close OUT;

  my $xsada=0
}
###############################################
sub splitToZFs{
  my ($oAln, $hRet) = @_;
  my ($zfStart,$current_pos,$zfNum,$uncertain,$some_ambiguity) = (0,1,1,0,0);

  for my $pc_id($oAln->consensus_conservation()){
      if ($pc_id < 80){
        $uncertain = 1;
        $some_ambiguity = 1;
      }
      ## Each ZF
      unless (++$current_pos % 84){
        $$hRet{$zfNum++} = $uncertain?"U"x84:substr($oAln->consensus_string(80),$zfStart,84);
        ($zfStart,$uncertain) = ($current_pos,0);
      }
  }

  return($some_ambiguity);
  #$$hRet{$zfNum++} = $uncertain?"U"x84:substr($oAln->consensus_string(80),$zfStart,84);
}

###############################################
sub makeFASTAs{
	my $stub = shift; 

	open  LHSFA, '>', $stub.".LHS.fa";
        print LHSFA ">lhs\n";
        print LHSFA "CACAGCCGTAATGACAAAACCAAAGGTCAAGAGATCAAAGAAAGGTCCAAACTCTTGAATAAAAGGACATGGCAGAGGGAGATT\n";
        print LHSFA "TCAAGGGCCTTTTCTAGCCCACCCAAAGGACAAATGGGGAGCTGTAGAGTGGGAAAAAGAATAATGGAAGAAGAGTCCAGAACA\n";
        print LHSFA "GGCCAGAAAGTGAATCCAGGGAACACAGGCAAATTATTTGTGGGGGTAGGAATCTCAAGAATTGCAAAAGTCAAGTATGGAGAG\n";
	close LHSFA;

	open  RHSFA, '>', $stub.".RHS.fa";
        print RHSFA ">rhs\n";
        print RHSFA "GATGAGTAAGTCATTAGTAATAAAACCTCATCTCAATAGCCACAAAAAGACAAATGTGGTCACCACACACTTGCACACCCCAGC\n";
        print RHSFA "TGTGAGGTGGCTTCAGCGGAAGTCTGCTGACCCCTTATATTCCCCGAGAGTATAAAGAGATCGGAAATAACTGATTAAACAAAT\n";
        print RHSFA "CCGCCACTTTCATGACTAGAGATGAGGAAGAACAAGGGATAGTTCTGTAAGTGTTCGGGGGACATCAGCATGTGTGGTTCTTTC\n";
	close RHSFA;

	open  ZFAFA, '>', $stub.".ZFA.fa";
        print ZFAFA ">zfa\n";
        print ZFAFA "TGTGGACAAGGTTTCAGTGTTAAATCAGATGTTATTACACACCAAAGGACACATACAGGGGAGAAGCTCTACGTCTGCAGGGAG\n";
	close ZFAFA;

	open  ZFFA, '>', $stub.".ZF.fa";
        print ZFFA ">zf\n";
        print ZFFA "TGTGGGCGGGGCTTTAGCCGGCAGTCAGTCCTCCTCACTCACCAGAGGAGACACACAGGGGAGAAGCCCTATGTCTGCAGGGAG\n";
	close ZFFA;

}

