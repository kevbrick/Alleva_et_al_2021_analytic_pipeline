#!/usr/bin/env perl
use strict;
use Getopt::Long;

GetOptions ('c=s'      => \(my $child_alleles),
            'p=s'      => \(my $parental_alleles),
            'f=s'      => \(my $prdm9_alleles_fasta),
            'o=s'      => \(my $output_file),
            'nohom+'   => \(my $do_not_test_hom),
            'nohet+'   => \(my $do_not_test_het),
            'nowarn+'  => \(my $no_warnings),
            'nohead+'  => \(my $no_header),
            'legdets+' => \(my $legacy_detail_names),
            't+'       => \(my $test_mode),
            'h+'       => \(my $help));

#-------------------------------------------------------------------------------------------------
showHelp () if ($help);

## TEST MODE OVERRIDES SETTINGS
if ($test_mode){
  open FA, '>', 'test.PRDM9TS.fa';
  makeTestFASTA();
  close FA;
  $prdm9_alleles_fasta = "test.PRDM9TS.fa";
  $output_file         = "test.PRDM9TS.out";
  $child_alleles       = "";
  $parental_alleles    = "";
  $do_not_test_hom     = "";
  $do_not_test_het     = "";
  $no_warnings++;
}

#-- CHECK ARGS ----------------------------------------------------------------------------------------
my $errMsg;

$errMsg .= "Allele file not specified (--f)\n"                       unless ($prdm9_alleles_fasta);
$errMsg .= "Allele file does not exist (--f $prdm9_alleles_fasta)\n" unless (-e $prdm9_alleles_fasta);

my (%seqs, @arr_child, @arr_parent);
 
read_in_alleles($prdm9_alleles_fasta,\%seqs,\$errMsg);

if ($child_alleles){
  @arr_child = are_alleles_OK($child_alleles,\$errMsg);
}else{
  print STDERR "** WARNING ** : Child allele(s) not provided (--c); Testing all alleles in FASTA file ($prdm9_alleles_fasta).\n" unless ($no_warnings);
  @arr_child = sort keys(%seqs);
}

my @arr_parent;
if ($parental_alleles){
  @arr_parent = are_alleles_OK($parental_alleles,\$errMsg);
}else{
  print STDERR "** WARNING ** : Parental allele(s) not provided (--p); Testing all alleles in FASTA file ($prdm9_alleles_fasta).\n" unless ($no_warnings);
  @arr_parent = sort keys(%seqs);
}

$errMsg .= "No child alleles found\n"  unless ($arr_child[0]);
$errMsg .= "No parent alleles found\n" unless ($arr_parent[0]);
showHelp ($errMsg) if ($errMsg);
#-- END CHECK ARGS ----------------------------------------------------------------------------------------

if ($output_file){
  open TSOUT, '>', $output_file;
}else{
  open TSOUT, '>&', STDOUT;
}

print_output_header() unless ($no_header);

for my $child(@arr_child){
  for my $i(0..$#arr_parent){
      for my $j($i..$#arr_parent){
        if ($arr_parent[$i] ne $child && $arr_parent[$j] ne $child){
          infer_template_switches($child, $arr_parent[$i], $arr_parent[$j]);
        }
    }
  }
}

close TSOUT; 

if ($test_mode){
  my $md5 = `md5sum $output_file`;
  $md5 =~ /^(\S+)\s+.+$/;
  if ($1 eq "9f47db08a7fc1e11309398830000f143"){
    print STDERR "TEST STATUS: OK\n"
  }else{
    print STDERR "TEST STATUS: MD5 does not match ... there might be something wrong. \n"
  }
}

################################################################################
sub showHelp {
  my ($msg) = @_;

  print "\n";
  print "infer_PRDM9_template_switches\n\n";
  print "------------------------------------------------------------------------------\n\n";
  print "Script to infer template switches that can give rise to a PRDM9 allele from\n";
  print "others (from: Alleva et al. 2021; under review). \n\n";
  
  if ($msg){  
    print "==============================================================================\n";
    print "ERRORS: \n";
    print "------------------------------------------------------------------------------\n";
    print $msg."\n";
  }
  
  print "- Arguments ------------------------------------------------------------------\n";
  print "--f       FASTA file (REQUIRED): Should only contain DNA sequences of PRDM9 \n";
  print "                                 ZFs for each allele. (REQUIRED)\n\n";
  print "--c       Child allele(s)      : Can be a single allele, comma-separated list,\n";
  print "                                 or a file with a list of IDs (one per line)\n";
  print "                                 *Use all sequences in FASTA if --c is null\n\n";
  print "--p       Parent allele(s)     : Can be a single allele, comma-separated list, \n";
  print "                                 or a file with a list of IDs (one per line)\n";
  print "                                 *Use all sequences in FASTA if --p is null\n\n";
  print "--o       Output file (default is STDOUT)\n\n";
  print "--nohom   Don't look for template switches if parental alleles are identical.\n\n";
  print "--nohet   Don't look for template switches if parental alleles differ.\n\n";
  print "--nohead  Don't print header on output file.\n\n";
  print "--t       TEST MODE - needs test.fa\n\n";
  print "--h       Show this help\n\n";
  print "------------------------------------------------------------------------------\n";
  print "contact : Kevin Brick (kevin.brick\@nih.gov)\n";
  print "------------------------------------------------------------------------------\n";
  
  exit;
}

################################################################################
sub read_in_alleles {
  my ($prdm9FA, $dets, $msg) = @_;
  
  my ($id, $seq);
  
  open IN, $prdm9FA;
  while (<IN>){
    chomp;
    if ($_ =~ /^\>(\S+)$/){
      $$dets{$id} = $seq if ($id);
      $id = $1;
      $seq = "";
    }else{
      if ($id){
        $seq .= $_;
      }
    }
    
    $$dets{$id} = $seq if ($seq);
  }
  $$msg .= "No sequences found. Is ".$prdm9FA." a FASTA file?\n" unless ($id);
}

################################################################################
sub are_alleles_OK {
  my ($alleles, $msg) = @_;
  
  my (@arr_alleles, @ret_alleles);
  
  if (-e $alleles){
    open ALLELES, $alleles; 
    while (<ALLELES>){
      chomp;
      push @arr_alleles, $_;
    }
    close ALLELES;
  }else{
    @arr_alleles = split(/,/,$alleles);
  }
  
  for my $allele (@arr_alleles){
    if ($seqs{$allele}){
        push @ret_alleles, $allele;
      }else{
        $$msg .= "Invalid allele [$allele]. NOT found in FASTA file.\n";
      }
  }
  
  return (@ret_alleles);
}

################################################################################
sub infer_template_switches{
  my ($id_child, $id_p1, $id_p2) = @_;

  my $recombStr;
  my $seq_child = $seqs{$id_child};
  
  if ($id_p1 eq $id_p2){
    unless ($do_not_test_hom){
      my $seq_parents = $seqs{$id_p1};
      $recombStr = rec($seq_parents,$seq_child,0)."\n";
      printTgtPatternReduced("HOM", $id_child, $id_p1, $id_p2, $recombStr,$seq_parents,$seq_child);
    }
  }else{
    unless ($do_not_test_het){
      my $seq_parents_12 = $seqs{$id_p1}."____".$seqs{$id_p2};
      my $seq_parents_21 = $seqs{$id_p2}."____".$seqs{$id_p1};
      $recombStr = rec($seq_parents_12,$seq_child,0)."\n";
      printTgtPatternReduced("HET", $id_child, $id_p1, $id_p2,$recombStr,$seq_parents_12,$seq_child);
      $recombStr = rec($seq_parents_21,$seq_child,0)."\n";
      printTgtPatternReduced("HET", $id_child, $id_p2, $id_p1,$recombStr,$seq_parents_21,$seq_child);
    }
  }
}

################################################################################
sub printTgtPatternReduced{
  my ($type, $name_child, $name_p1, $name_p2, $str, $seq, $qry) = @_;

  my $prevX = '';
  my $segStr = $str; chomp $segStr;

  my @seqs = split(/_+/,$seq);
  $seq =~ s/_//g;

  my $allele1Len = floor(length(@seqs[0])/84);
  my $allele2Len = floor(length(@seqs[1])/84);
  my $queryLen   = floor(length($qry)/84);

  my %schema;
  my @schema  = split("",("0" x $queryLen));

  my $prev=0;
  my $parental_allele = 1;
  while ($str =~ s/^.+?(\d+)\-(\d+)//){
    my ($ffrom,$fto) = (($1+1)/84,($2+1)/84);
    my ($from,$to) = (floor($ffrom),floor($fto));

    #my $parental_allele = ($to <= $allele1Len)?1:2;

    if ($ffrom == int($ffrom)){
      $prev += 1;
    }

    for my $i($prev..$prev+($to-$from)){
      $schema[$i] = $schema[$i]?3:$parental_allele;
    }
    $prev=$prev+($to-$from);
    
    if ($name_p1 ne $name_p2){
      $parental_allele = $parental_allele == 1?2:1;
    }
  }
  
  my $recType;
  if ($legacy_detail_names){
    $recType = ($name_p1 eq $name_p2)?"Monoallelic":"Recombinant";
  }else{
    $recType = ($name_p1 eq $name_p2)?"PARENT_HOM":"PARENT_HET";
  }
  
  my @nShifts = split(/::/,$segStr); 
  if ($segStr eq "NA"){
    print TSOUT join("\t",$name_child,$name_p1,$name_p2,$recType,"NA","NA","NA")."\n";
  }else{
    print TSOUT join("\t",$name_child,$name_p1,$name_p2,$recType,$segStr,join(":",@schema),$#nShifts)."\n";
  }
}  

################################################################################
sub pad{
  my ($toPad,$n) = @_;
  return($toPad.(" " x ($n - length($toPad))));
}

################################################################################
sub rec{
  my ($p1,$res,$start) = @_;

  my $lhsLen = lhsMatch($p1,$res);
  return "NA" unless ($lhsLen > -1);

  my $recombDets = ":".$start."-".($start+$lhsLen-1).":";

  if ($lhsLen == length($res)){
    return $recombDets.":END";
  }

  ## Now, ask if the remainder can be found in the other allele?
  my $remaining_query_seq = substr($res,$lhsLen,length($res));

  while ($remaining_query_seq){
    if ($p1 =~ /^(.*)($remaining_query_seq)/){
      #THIS IS A RECOMBINANT ALLELE
      return $recombDets.":".($start+length($1))."-".($start+length($1)+length($remaining_query_seq)-1).":END";
    }else{
      my ($tPos,$qPos,$matchLen);
      ($tPos,$qPos,$matchLen,$remaining_query_seq) = nextLHSmatch($p1,$remaining_query_seq);
      #($tPos,$qPos,$matchLen,$remaining_query_seq) = getNextLHSmatch($p1,$remaining_query_seq);

      return "NA" unless ($tPos > -1);

      $recombDets .= ":".$tPos."-".($tPos+$matchLen-1).":";
    }
  }
}

################################################################################
sub lhsMatch {
  my ($tgt,$query) = @_;
  ## Find longest LHS match using binary matching
  $_ = $tgt ^ $query;
  s/./ord $& ? "^" : " "/ge;

  $_ =~ /^([^\^]+?)\^/;

  #print STDERR "L:".substr($query,0,length($1))."\n";
  #return(length($1)-1) if ($1);
  return(length($1)) if ($1);
  return -1;
}

################################################################################
sub floor{
  my $num = shift;
  return(int(sprintf("%4i",$num)));
}

################################################################################
sub round{
  my $num = shift;
  return(int(sprintf("%4i",$num+0.5)));
}

################################################################################
sub nextLHSmatch{
  my ($target_seq,$query_seq) = @_;

  my ($subQ,$iteration);
  my $cc = 999;

  while ($cc > 1){
    $cc = 0;
    $subQ = substr($query_seq,0,9+($iteration++));
    ++$cc while($target_seq =~ /($subQ)/g);
  }

  ## if there are multiple matches to a short segment, use the first
  if ($cc == 0){
    $subQ = substr($query_seq,0,9+($iteration-2));
  }

  $target_seq =~ s/^(.*?)($subQ)/$2/g;

  my $tgtPos = length($1)+1;

  my $lhsMaxLen = lhsMatch($target_seq,$query_seq);

  return (-1,-1,-1,'') if ($lhsMaxLen <= 10);

  my $remaining_seq = substr($query_seq, $lhsMaxLen, length($query_seq) - $lhsMaxLen);

  return($tgtPos,0,$lhsMaxLen,$remaining_seq);

}

################################################################################
sub makeTestFASTA{
  my %prdm9;
  
  $prdm9{A} = "TGTGGACAAGGTTTCAGTGTTAAATCAGATGTTATTACACACCAAAGGACACATACAGGGGAGAAGCTCTACGTCTGCAGGGAGTGTGGGCGGGGCTTTAGCTGGAAGTCACACCTCCTCATTCACCAGAGGATACACACAGGGGAGAAGCCCTATGTCTGCAGGGAGTGTGGGCGGGGCTTTAGCTGGCAGTCAGTCCTCCTCACTCACCAGAGGACACACACAGGGGAGAAGCCCTATGTCTGCAGGGAGTGTGGGCGGGGCTTTAGCCGGCAGTCAGTCCTCCTCACTCACCAGAGGAGACACACAGGGGAGAAGCCCTATGTCTGCAGGGAGTGTGGGCGGGGCTTTAGCCGGCAGTCAGTCCTCCTCACTCACCAGAGGAGACACACAGGGGAGAAGCCCTATGTCTGCAGGGAGTGTGGGCGGGGCTTTAGCTGGCAGTCAGTCCTCCTCAGTCACCAGAGGACACACACAGGGGAGAAGCCCTATGTCTGCAGGGAGTGTGGGCGGGGCTTTAGCTGGCAGTCAGTCCTCCTCACTCACCAGAGGACACACACAGGGGAGAAGCCCTATGTCTGCAGGGAGTGTGGGCGGGGCTTTAGCAATAAGTCACACCTCCTCAGACACCAGAGGACACACACAGGGGAGAAGCCCTATGTCTGCAGGGAGTGTGGGCGGGGCTTTCGCGATAAGTCACACCTCCTCAGACACCAGAGGACACACACAGGGGAGAAGCCCTATGTCTGCAGGGAGTGTGGGCGGGGCTTTAGAGATAAGTCAAACCTCCTCAGTCACCAGAGGACACACACAGGGGAGAAGCCCTATGTCTGCAGGGAGTGTGGGCGGGGCTTTAGCAATAAGTCACACCTCCTCAGACACCAGAGGACACACACAGGGGAGAAGCCCTATGTCTGCAGGGAGTGTGGGCGGGGCTTTCGCAATAAGTCACACCTCCTCAGACACCAGAGGACACACACAGGGGAGAAGCCCTACGTCTGCAGGGAGTGTGGGCGGGGCTTTAGCGATAGGTCAAGCCTCTGCTATCACCAGAGGACACACACAGGGGAGAAGCCCTACGTCTGCAGGGAG";

  $prdm9{B} = "TGTGGACAAGGTTTCAGTGTTAAATCAGATGTTATTACACACCAAAGGACACATACAGGGGAGAAGCTCTACGTCTGCAGGGAGTGTGGGCGGGGCTTTAGCTGGAAGTCACACCTCCTCATTCACCAGAGGATACACACAGGGGAGAAGCCCTATGTCTGCAGGGAGTGTGGGCGGGGCTTTAGCTGGCAGTCAGTCCTCCTCACTCACCAGAGGACACACACAGGGGAGAAGCCCTATGTCTGCAGGGAGTGTGGGCGGGGCTTTAGCCGGCAGTCAGTCCTCCTCACTCACCAGAGGAGACACACAGGGGAGAAGCCCTATGTCTGCAGGGAGTGTGGGCGGGGCTTTAGCCGGCAGTCAGTCCTCCTCACTCACCAGAGGAGACACACAGGGGAGAAGCCCTATGTCTGCAGGGAGTGTGGGCGGGGCTTTAGCTGGCAGTCAGTCCTCCTCACTCACCAGAGGACACACACAGGGGAGAAGCCCTATGTCTGCAGGGAGTGTGGGCGGGGCTTTAGCTGGCAGTCAGTCCTCCTCACTCACCAGAGGACACACACAGGGGAGAAGCCCTATGTCTGCAGGGAGTGTGGGCGGGGCTTTAGCAATAAGTCACACCTCCTCAGACACCAGAGGACACACACAGGGGAGAAGCCCTATGTCTGCAGGGAGTGTGGGCGGGGCTTTCGCGATAAGTCACACCTCCTCAGACACCAGAGGACACACACAGGGGAGAAGCCCTATGTCTGCAGGGAGTGTGGGCGGGGCTTTAGAGATAAGTCAAACCTCCTCAGTCACCAGAGGACACACACAGGGGAGAAGCCCTATGTCTGCAGGGAGTGTGGGCGGGGCTTTAGCAATAAGTCACACCTCCTCAGACACCAGAGGACACACACAGGGGAGAAGCCCTATGTCTGCAGGGAGTGTGGGCGGGGCTTTCGCAATAAGTCACACCTCCTCAGACACCAGAGGACACACACAGGGGAGAAGCCCTACGTCTGCAGGGAGTGTGGGCGGGGCTTTAGCGATAGGTCAAGCCTCTGCTATCACCAGAGGACACACACAGGGGAGAAGCCCTACGTCTGCAGGGAG";
  
  $prdm9{N} = "TGTGGACAAGGTTTCAGTGTTAAATCAGATGTTATTACACACCAAAGGACACATACAGGGGAGAAGCTCTACGTCTGCAGGGAGTGTGGGCGGGGCTTTAGCTGGAAGTCACACCTCCTCATTCACCAGAGGATACACACAGGGGAGAAGCCCTATGTCTGCAGGGAGTGTGGGCGGGGCTTTAGCTGGCAGTCAGTCCTCCTCACTCACCAGAGGACACACACAGGGGAGAAGCCCTATGTCTGCAGGGAGTGTGGGCGGGGCTTTAGCCGGCAGTCAGTCCTCCTCACTCACCAGAGGAGACACACAGGGGAGAAGCCCTATGTCTGCAGGGAGTGTGGGCGGGGCTTTAGCCGGCAGTCAGTCCTCCTCACTCACCAGAGGAGACACACAGGGGAGAAGCCCTATGTCTGCAGGGAGTGTGGGCGGGGCTTTAGCTGGCAGTCAGTCCTCCTCAGTCACCAGAGGACACACACAGGGGAGAAGCCCTATGTCTGCAGGGAGTGTGGGCGGGGCTTTAGCTGGCAGTCAGTCCTCCTCACTCACCAGAGGACACACACAGGGGAGAAGCCCTATGTCTGCAGGGAGTGTGGGCGGGGCTTTAGCAATAAGTCACACCTCCTCAGACACCAGAGGACACACACAGGGGAGAAGCCCTATGTCTGCAGGGAGTGTGGGCGGGGCTTTCGCGATAAGTCACACCTCCTCAGACACCAGAGGACACACACAGGGGAGAAGCCCTATGTCTGCAGGGAGTGTGGGCGGGGCTTTAGAGATAAGTCAAACCTCCTCAGTCACCAGAGGACACACACAGGGGAGAAGCCCTATGTCTGCAGGGAGTGTGGGCGGGGCTTTAGCAATAAGTCACACCTCCTCAGACACCAGAGGACACACACAGGGGAGAAGCCCTACGTCTGCAGGGAGTGTGGGCGGGGCTTTAGCGATAGGTCAAGCCTCTGCTATCACCAGAGGACACACACAGGGGAGAAGCCCTACGTCTGCAGGGAG";
  
  $prdm9{C} = "TGTGGACAAGGTTTCAGTGTTAAATCAGATGTTATTACACACCAAAGGACACATACAGGGGAGAAGCTCTACGTCTGCAGGGAGTGTGGGCGGGGCTTTAGCTGGAAGTCACACCTCCTCATTCACCAGAGGATACACACAGGGGAGAAGCCCTATGTCTGCAGGGAGTGTGGGCGGGGCTTTAGCTGGCAGTCAGTCCTCCTCACTCACCAGAGGACACACACAGGGGAGAAGCCCTATGTCTGCAGGGAGTGTGGGCGGGGCTTTAGCCGGCAGTCAGTCCTCCTCACTCACCAGAGGAGACACACAGGGGAGAAGCCCTATGTCTGCAGGGAGTGTGGGCGGGGCTTTAGCCGGCAGTCAGTCCTCCTCACTCACCAGAGGAGACACACAGGGGAGAAGCCCTATGTCTGCAGGGAGTGTGGGCGGGGCTTTAGCTGGCAGTCAGTCCTCCTCACTCACCAGAGGACACACACAGGGGAGAAGCCCTATGTCTGCAGGGAGTGTGGGCGGGGCTTTAGCTGGCAGTCAGTCCTCCTCACTCACCAGAGGACACACACAGGGGAGAAGCCCTATGTCTGCAGGGAGTGTGGGCGGGGCTTTAGCAATAAGTCACACCTCCTCAGACACCAGAGGACACACACAGGGGAGAAGCCCTATGTCTGCAGGGAGTGTGGGCGGGGCTTTAGAGATAAGTCACACCTCCTCAGTCACCAGAGGACACACACAGGGGAGAAGCCCTATGTCTGCAGGGAGTGTGGGCGGGGCTTTAGAGATAAGTCAAACCTCCTCAGTCACCAGAGGACACACACAGGGGAGAAGCCCTATGTCTGCAGGGAGTGTGGGCGGGGCTTTAGCTGGCAGTCAGTCCTCCTCAGACACCAGAGGACACACACAGGGGAGAAGCCCTATGTCTGCAGGGAGTGTGGGCGGGGCTTTAGAGATAAGTCAAACCTCCTCAGTCACCAGAGGACACACACAGGGGAGAAGCCCTATGTCTGCAGGGAGTGTGGGCGGGGCTTTCGCAATAAGTCACACCTCCTCAGACACCAGAGGACACACACAGGGGAGAAGCCCTACGTCTGCAGGGAGTGTGGGCGGGGCTTTAGCGATAGGTCAAGCCTCTGCTATCACCAGAGGACACACACAGGGGAGAAGCCCTACGTCTGCAGGGAG";
  
  $prdm9{L4} = "TGTGGACAAGGTTTCAGTGTTAAATCAGATGTTATTACACACCAAAGGACACATACAGGGGAGAAGCTCTACGTCTGCAGGGAGTGTGGGCGGGGCTTTAGCTGGAAGTCACACCTCCTCATTCACCAGAGGATACACACAGGGGAGAAGCCCTATGTCTGCAGGGAGTGTGGGCGGGGCTTTAGCTGGCAGTCAGTCCTCCTCACTCACCAGAGGACACACACAGGGGAGAAGCCCTATGTCTGCAGGGAGTGTGGGCGGGGCTTTAGCCGGCAGTCAGTCCTCCTCACTCACCAGAGGAGACACACAGGGGAGAAGCCCTATGTCTGCAGGGAGTGTGGGCGGGGCTTTAGCCGGCAGTCAGTCCTCCTCACTCACCAGAGGAGACACACAGGGGAGAAGCCCTATGTCTGCAGGGAGTGTGGGCGGGGCTTTAGCTGGCAGTCAGTCCTCCTCACTCACCAGAGGACACACACAGGGGAGAAGCCCTATGTCTGCAGGGAGTGTGGGCGGGGCTTTAGCTGGCAGTCAGTCCTCCTCACTCACCAGAGGACACACACAGGGGAGAAGCCCTATGTCTGCAGGGAGTGTGGGCGGGGCTTTAGCTGGCAGTCAGTCCTCCTCACTCACCAGAGGACACACACAGGGGAGAAGCCCTATGTCTGCAGGGAGTGTGGGCGGGGCTTTAGCCGGCAGTCAGTCCTCCTCACTCACCAGAGGAGACACACAGGGGAGAAGCCCTATGTCTGCAGGGAGTGTGGGCGGGGCTTTAGCCGGCAGTCAGTCCTCCTCACTCACCAGAGGAGACACACAGGGGAGAAGCCCTATGTCTGCAGGGAGTGTGGGCGGGGCTTTAGCTGGCAGTCAGTCCTCCTCACTCACCAGAGGACACACACAGGGGAGAAGCCCTATGTCTGCAGGGAGTGTGGGCGGGGCTTTAGCAATAAGTCACACCTCCTCAGACACCAGAGGACACACACAGGGGAGAAGCCCTATGTCTGCAGGGAGTGTGGGCGGGGCTTTAGAGATAAGTCACACCTCCTCAGTCACCAGAGGACACACACAGGGGAGAAGCCCTATGTCTGCAGGGAGTGTGGGCGGGGCTTTAGAGATAAGTCAAACCTCCTCAGTCACCAGAGGACACACACAGGGGAGAAGCCCTATGTCTGCAGGGAGTGTGGGCGGGGCTTTAGCTGGCAGTCAGTCCTCCTCAGACACCAGAGGACACACACAGGGGAGAAGCCCTATGTCTGCAGGGAGTGTGGGCGGGGCTTTAGAGATAAGTCAAACCTCCTCAGTCACCAGAGGACACACACAGGGGAGAAGCCCTATGTCTGCAGGGAGTGTGGGCGGGGCTTTCGCAATAAGTCACACCTCCTCAGACACCAGAGGACACACACAGGGGAGAAGCCCTACGTCTGCAGGGAGTGTGGGCGGGGCTTTAGCGATAGGTCAAGCCTCTGCTATCACCAGAGGACACACACAGGGGAGAAGCCCTACGTCTGCAGGGAG";
  
  for my $i (sort keys(%prdm9)) {
    print FA ">$i\n$prdm9{$i}\n";
  }
  close FA;

}

################################################################################
sub print_output_header {
  print TSOUT "## Output of template switching (TS) algorithm\n";
  print TSOUT "## -------------------------------------------------------------------------------\n";
  print TSOUT "## Column details: \n";
  print TSOUT "## -------------------------------------------------------------------------------\n";
  print TSOUT "## child_allele               Child PRDM9 allele\n";
  print TSOUT "##                            \n";
  print TSOUT "## starting_parental_allele   Parental PRDM9 allele that serves as the starting\n";
  print TSOUT "##                            allele for the TS algorithm \n";
  print TSOUT "##                            \n";
  print TSOUT "## other_parental_allele      The second PRDM9 parental allele. This can be the \n";
  print TSOUT "##                            same as the first allele. \n";
  print TSOUT "##                            \n";
  print TSOUT "## ts_type                    HOM = switches are between identical alleles\n";
  print TSOUT "##                            HET = switches are between different alleles\n";
  print TSOUT "##                            \n";
  print TSOUT "## ts_details                 TS path details. Each homology tract is described as\n";
  print TSOUT "##                            :{start}-{end}:\n";
  print TSOUT "##                            The end of a tract indicates a template switch \n";
  print TSOUT "##                            (unless it is the final tract, shown as END). \n";
  print TSOUT "##                            Example:\n";
  print TSOUT "##                            Human PRDM9-N can arise from one TS between PRDM9-A\n";
  print TSOUT "##                            and PRDM9-C. This event is described as:\n";
  print TSOUT "##                            :0-910::2175-2271:END	\n";
  print TSOUT "##                            Tract 1 = 0-910 from PRDM9-A\n";
  print TSOUT "##                            -- Template switch -- \n";
  print TSOUT "##                            Tract 2 = 2175-2271 from PRDM9-C\n";
  print TSOUT "##                            \n";
  print TSOUT "## parentage_map              Schematic of the parental origins of the child allele. \n";
  print TSOUT "##                              1 = derived from parent 1\n";
  print TSOUT "##                              2 = derived from parent 2\n";
  print TSOUT "##                              3 = location of template switch\n";
  print TSOUT "##                                  NOTE: Resolution is 1 ZF; if >1 TS occurs in a \n";
  print TSOUT "##                                        single ZF, it is not shown. This can result\n";
  print TSOUT "##                                        the same parental origin both sides of an \n";
  print TSOUT "##                                        apparent TF. This can be checked in ts_details.\n";
  print TSOUT "##                             E.G. TS events from PRDM9-A (1) and PRDM9-C (2) to make PRDM9-N\n";
  print TSOUT "##                                  1:1:1:1:1:1:1:1:1:1:3:2:2\n";
  print TSOUT "##                                  |From Pr-A        | ^ |From Pr-C\n";
  print TSOUT "##                                                      |\n";
  print TSOUT "##                                                      |\n";
  print TSOUT "##                                                     T.S\n";
  print TSOUT "##                            \n";
  print TSOUT "## num_switches               Number of template switches required \n";
  print TSOUT "##                          \n";
  print TSOUT "## NOTE: If no TS can yield the child allele, the last three columns will be \"NA\"\n";  
    
  print TSOUT join("\t","child_allele",
                       "starting_parental_allele",
                       "other_parental_allele",
                       "ts_type",
                       "ts_details",
                       "ts_simplified",
                       "num_switches")."\n";
}