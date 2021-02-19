use strict;
use String::LCSS;
use Getopt::Long;

GetOptions ('a=s' => \(my $allelesFile = 'PrZFA_allele_codes.txt'),
            'z=s' => \(my $zfFile      = 'PrZFA_ZF_sequences.txt'),
            'c=s' => \(my $allelesToCheck),
            'p=s' => \(my $parentalAllele),
            'o=s' => \(my $out         = 'PrZFA_recombinants.txt'));

## If alleles are in a file
my @selectAlleles;

if ($allelesToCheck){
  if (-e $allelesToCheck){
    open ALLELE, $allelesToCheck;
    while(<ALLELE>){chomp; push @selectAlleles, $_};
    close ALLELE;
  }else{
    @selectAlleles = ($allelesToCheck);
  }
}

## Read alleles and order them
my (%zfs, %alleles, %recomb);

readHumanZFs($zfFile, \%zfs);
readHumanAlleles($allelesFile, \%zfs, \%alleles);

$alleles{'N'} = $alleles{'Av:s:0053:M1S:A-A'} if($alleles{'Av:s:0053:M1S:A-A'});

my (@az, @L, @M, @B, @v, @Av, @other);
for my $allele(sort keys %alleles){
  if ($allele =~ /^[A-Z]$/i)  {push @az, $allele ; next};
  if ($allele =~ /^pt[A-Z]$/i){push @az, $allele ; next};
  if ($allele =~ /^L\d+$/)    {push @L,  $allele ; next};
  if ($allele =~ /^M\d+$/)    {push @M,  $allele ; next};
  if ($allele =~ /^baudat/)   {push @B,  $allele ; next};

  if ($allele =~ /^Av:/)      {push @Av,  $allele };
  if ($allele =~ /v:/)        {push @v,  $allele ; next};

  if ($allele =~ /(MOL|CST|DOM|SPR|UNK)_(\S+)/) {push @M, $allele; next};
  push @other, $allele;
}

#Set arrays of parental alleles and alleles to check
my @checkAlleles = $allelesToCheck?@selectAlleles:(@az, @L, @M, @B);
my @allAlleles   = $parentalAllele?($parentalAllele):(@az, @L, @M, @B, @v, @other);

open OUT, '>', $out;
print OUT    join("\t","PrZFA","Parental","event","nsteps")."\n";
print STDERR join("\t","PrZFA","Parental","event","nsteps")."\n";

for my $queryAllele(@checkAlleles){
  for my $parentAllele(@allAlleles){

    next if ($parentAllele eq $queryAllele);

    $recomb{$queryAllele}->{$parentAllele} = rec($alleles{$parentAllele}->{'seq'},$alleles{$queryAllele}->{'seq'}, 0);
    my $nRecs = ($recomb{$queryAllele}->{$parentAllele} =~ /NOREC/)?0:(($recomb{$queryAllele}->{$parentAllele} =~ tr/\://) - 2)/2;
    print OUT    join("\t",$queryAllele,$parentAllele,$recomb{$queryAllele}->{$parentAllele},$nRecs>0?$nRecs:0)."\n";
    print STDERR join("\t",$queryAllele,$parentAllele,$recomb{$queryAllele}->{$parentAllele},$nRecs>0?$nRecs:0)."\n";

  }
}
close OUT;

################################################################################
sub rec{
  my ($p1,$res,$start) = @_;

  my $lhsLen = lhsMatch($p1,$res);
  return "NORECOMBINANT" unless ($lhsLen>0);

  ## Now, ask if the remainder can be found in the other allele?
  my $rest = substr($res,$lhsLen,length($res));

  my $recombDets = ":".$start."-".($start+$lhsLen).":";

  if ($p1 =~ /^(.*)($rest)/){
    #THIS IS A RECOMBINANT ALLELE
    return $recombDets.":".($start+length($1))."-".($start+length($1)+length($rest)).":END";
  }else{
    my $nextPos = getNextLHSmatch($p1,$rest);

    return "NORECOMBINANT" unless ($nextPos);

    my $p1New = substr($p1,$nextPos,length($p1)-$nextPos);
    return $recombDets.rec($p1New,$rest,$nextPos);

  }
}

################################################################################
sub lhsMatch {
  my ($tgt,$query) = @_;
  ## Find longest LHS match using binary matching
  $_ = $tgt ^ $query;
  s/./ord $& ? "^" : " "/ge;

  $_ =~ /^(.+?)\^/;

  return(length($1)-1) if ($1);
  return 0;
}

################################################################################
sub getNextLHSmatch {
  my ($tgt, $query) = @_;

  my $queryPos = -1;
  my (@result,$len,$p1New, $tgtPos);

  my $lcssQ = $query;

  while ($queryPos != 0){

    ## Pain in the arse ... index position depends on which is the longer sequence !!!
    my ($lcss, $idxShorter, $idxLonger) = String::LCSS::lcss ( $tgt , $lcssQ );

    if (length($tgt) < length($lcssQ)){
      $queryPos = $idxLonger;
      $tgtPos   = $idxShorter;
    }else{
      $queryPos = $idxShorter;
      $tgtPos   = $idxLonger;
    }

    $len = length($lcss)-1;

    if ($queryPos){
      $lcssQ = substr($lcssQ,0,$queryPos)
    }else{
      return 0 if ($len <= 3);
      return($tgtPos);
    }
  }
  return 0;
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
