use strict;
use String::LCSS_XS;
use Getopt::Long;

GetOptions ('n=s'    => \(my $names),
            'mode=s' => \(my $mode = 'test'),
            'sP1=s'  => \(my $pSeq1),
            'sP2=s'  => \(my $pSeq2),
            'sC=s'   => \(my $seqKid),
            'all=s'  => \(my $prdm9_alleles));

my (%seqs, %codes);

if ($prdm9_alleles && -e $prdm9_alleles){
  open PIN, $prdm9_alleles;
  while (<PIN>){
    chomp;
    next if ($_ =~ /^(#|ID)/);
    my ($ID,$sID,$zf_code,$pub,$in_pop,$dna_sequence,$aa,$dna_contact_aas,$ABDdist,$CBDdist,$ACtype) =split("\t",$_);
    $seqs{$sID} = $dna_sequence;
    $seqs{$ID} = $dna_sequence;
    
    $codes{$sID} = $zf_code;
    $codes{$ID}  = $zf_code;
  }
}

### TEST MODE :
### Rv:0061 from A & L20
if ($mode == 'test'){ 
  $names = "A,L20,Rv:0061";
  $mode = 3;
}

my ($nameP1, $nameP2, $nameChild) = split(/,/,$names);

$pSeq1  = $seqs{$nameP1} if ($seqs{$nameP1});
$pSeq2  = $seqs{$nameP2} if ($seqs{$nameP2});
$seqKid = $seqs{$nameChild} if ($seqs{$nameChild});

my $pSeq12 = $pSeq1."____".$pSeq2;
my $pSeq21 = $pSeq2."____".$pSeq1;

if ($mode ne 1 && $nameP1 eq $nameP2){
  print STDERR "********************************************************** \n";
  print STDERR "$nameP1 == $nameP2\n";
  print STDERR "MODE set to 1 (from $mode); Identical parental alleles ... \n";
  print STDERR "********************************************************** \n";
  $mode = 1;
}

die ("ERROR: No sequence for $nameP1\n")    unless($pSeq1);
die ("ERROR: No sequence for $nameP2\n")    unless($pSeq2);
die ("ERROR: No sequence for $nameChild\n") unless($seqKid);

my ($recombStr); 

if ($mode == 1 || $mode == 3){ ## Check ONLY PARENTAL ALLELES - NO RECOMBINATION
  if ($nameP1 eq $nameP2){
    checkForRecombinants($pSeq1,$seqKid,join("|p12p|",$nameP1,$nameP1));
    exit
  }else{
    checkForRecombinants($pSeq1,$seqKid,join("|p12p|",$nameP1,$nameP1));
  }
    
  checkForRecombinants($pSeq2,$seqKid,join("|p12p|",$nameP2,$nameP2));
  exit unless ($mode eq 3);
}

if ($mode == 2 || $mode == 3){ ## Check ONLY PARENTAL ALLELES - NO RECOMBINATION
  checkForRecombinants($pSeq12,$seqKid,join("|p12p|",$nameP1,$nameP2));
  checkForRecombinants($pSeq21,$seqKid,join("|p12p|",$nameP2,$nameP1));
  exit;
}

################################################################################
sub printQryPatternReduced{
  my ($type, $str, $seq) = @_;

  my $prevX = '';
  #print "**** TEMPLATE SWITCHES ****\n";
  
  my @seqs = split(/_+/,$seq);
  $seq =~ s/_//g;
    
  my $allele1Len = floor(length(@seqs[0])/84);
  my $allele2Len = floor(length(@seqs[1])/84);
  
  my %schema;
  my @schema1 = split("",("." x $allele1Len));
  my @schema2 = split("",("." x $allele2Len));
  $schema{1} = \@schema1;
  $schema{2} = \@schema2;
  
  #print join("\t",pad("from",4),"-",pad("to",4))."\n";
  my $cnt=0;
  while ($str =~ s/^.+?(\d+)\-(\d+)//){
    
    $cnt++;
    my ($from,$to) = (floor($1/84),floor($2/84));
    if ($to <= $allele1Len){
      my $s = $schema{1};
      @$s[$from..$to] = split("","1" x ($to-$from+1));
      #print join(" ",@$s)."\n";
      print "TS\t$type\t$cnt\t".join("\t",split(",",$names))."\t".join(":",@$s)."\n";
    }else{
      my $s = $schema{2};
      $from -= (1+$allele1Len);
      $to   -= (1+$allele1Len);
      @$s[$from..$to] = split("","2" x ($to-$from+1));
      #print join(" ",@$s)."\n";
      print "TS\t$type\t$cnt\t".join("\t",split(",",$names))."\t".join(":",@$s)."\n";
    }
  }
}  

################################################################################
sub printTgtPatternReduced{
  my ($type, $str, $seq, $qry) = @_;

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
  while ($str =~ s/^.+?(\d+)\-(\d+)//){
    my ($ffrom,$fto) = (($1+1)/84,($2+1)/84);
    my ($from,$to) = (floor($ffrom),floor($fto));

    my $parental_allele = ($to <= $allele1Len)?1:2;
    if ($ffrom == int($ffrom)){
      $prev += 1;
    }
    
    for my $i($prev..$prev+($to-$from)){
      $schema[$i] = $schema[$i]?3:$parental_allele;
    }
    $prev=$prev+($to-$from);
  }
  my ($name1,$name2) = split(/\|p12p\|/,$type);
  my $recType = ($name1 eq $name2)?"Monoallelic":"Recombinant";
  my @nShifts = split(/::/,$segStr); 
  print join("\t",$nameChild,$name1,$name2,$recType,$segStr,join(":",@schema),$#nShifts)."\n";
}  

################################################################################
sub pad{
  my ($toPad,$n) = @_;
  return($toPad.(" " x ($n - length($toPad))));
}

################################################################################
sub checkForRecombinants{
  my ($seq_parents, $seq_child, $type_name) = @_;
  
  $recombStr = rec($seq_parents,$seq_child,0)."\n";
      
  #printQryPatternReduced($type_name,$recombStr,$seq_parents);
  printTgtPatternReduced($type_name,$recombStr,$seq_parents,$seq_child);
}

################################################################################
sub rec{
  my ($p1,$res,$start) = @_;
  
  my $lhsLen = lhsMatch($p1,$res);
  return "NORECOMBINANT" unless ($lhsLen > -1);

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
      ($tPos,$qPos,$matchLen,$remaining_query_seq) = getNextLHSmatch($p1,$remaining_query_seq);

      return "NORECOMBINANT" unless ($tPos > -1);
      
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
sub getNextLHSmatch {
  my ($tgt, $query) = @_;

  my $queryPos = -1;
  my (@result,$len,$p1New, $tgtPos);

  my $lcssQ = $query;

  while ($queryPos != 0){

    # ## Pain in the arse ... index position depends on which is the longer sequence !!!
    # my ($lcss, $idxShorter, $idxLonger) = String::LCSS_XS::lcss ( $tgt , $lcssQ );
    # 
    # if (length($tgt) < length($lcssQ)){
    #   $queryPos = $idxLonger;
    #   $tgtPos   = $idxShorter;
    # }else{
    #   $queryPos = $idxShorter;
    #   $tgtPos   = $idxLonger;
    # }
    
    ## WAAAAAAAAAAAAAAAAAAAAAAAAAAAY faster !!!!!
    my ($lcss, $tgtPos, $queryPos) = String::LCSS_XS::lcss ( $tgt , $lcssQ );

    $len = length($lcss)-1;

    if ($queryPos){
      $lcssQ = substr($lcssQ,0,$queryPos)
    }else{
      my $lhsMaxLen = lhsMatch(substr($tgt,$tgtPos,length($tgt)-$tgtPos),$query);
      return (-1,-1,-1,'') if ($lhsMaxLen <= 10);

      my $remaining_seq = substr($query, $lhsMaxLen, length($query) - $lhsMaxLen);
      return($tgtPos,$queryPos,$lhsMaxLen,$remaining_seq);
    }
  }
  return (-1,-1,-1,'');
}

################################################################################
sub readHumanZFs{
  my ($zfFile, $zfRef) = @_;
  open IN, $zfFile;
  while (<IN>){
    chomp;
    my ($name,$oldcode,$src,$seq) = split(/\t/,$_);
    $$zfRef{$name} = $seq;
    $$zfRef{$seq} = $name;
  }
  close IN;
}

################################################################################
sub readHumanAlleles{
  my ($alleleFile, $zfRef, $alleleRef) = @_;
  open IN, $alleleFile;
  while (<IN>){
    chomp;
    my $seq;
    my ($name,$code) = split(/\t/,$_);

    my $c2 = $code;

    while ($code =~ /(.{2})/g){
      $seq .= $$zfRef{$1};
    }
    $$alleleRef{$name}->{code} = $code;
    $$alleleRef{$name}->{seq}  = $seq;
  }
  close IN;
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