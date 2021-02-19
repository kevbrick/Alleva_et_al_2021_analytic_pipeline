use strict;
use List::MoreUtils qw(uniq);

my %zfs;
my %alleles;
my %mouseAlleles;
my %pubAlleles;
my %mouseZFs;
my %done;
my %seqDone;

my $zfN   = 64;
my $nameN = 64;
my @seps  = ("_",";",".","-");
my $sep   = 0;

#####################################################
open PUB, "mouse_PrZFA_alleles_in_Pubs.txt";
while (<PUB>){
  chomp;
  next if ($_ =~ /name/);
  my ($id,$name) = split(/\t/,$_);
  $pubAlleles{$id} = $id."::".$name;
}
close PUB;
#####################################################

open my $IN, "-|", "sort -k1,1 mousePrZFA.tab";

while (<$IN>){
  chomp;
  my ($name,$seq) = split(/\t/,$_);
  #my $alleleName = $done{$name}++?"mm".chr($nameN):"mm".chr(++$nameN);
  my $alleleName = $pubAlleles{$name}?$pubAlleles{$name}:$name;

  ## Skip non alphanumerics
  $nameN = $nameN == 90?97:$nameN;
  if ($nameN == 122){
    $nameN = 64;
  }

  ## Skip non alphanumerics
  $zfN = $zfN == 90?97:$zfN;
  if ($zfN == 122){
    $zfN = 64;
    $sep++;
  }

  my $zf = $zfs{$seq}?$zfs{$seq}:$seps[$sep].chr(++$zfN);
  $zfs{$seq} = $zf;
  $mouseZFs{$zf} = $seq;

  $alleles{$alleleName}->{code} .= $zf;
  $alleles{$alleleName}->{seq}  .= $seq;
  $mouseAlleles{$alleleName}->{code} .= $zf;
  $mouseAlleles{$alleleName}->{seq}  .= $seq;
}
close $IN;

open ZFS, '>', 'PrZFA_ZF_sequences.mouse.txt';
for my $zf(uniq(sort keys(%mouseZFs))) {
  print ZFS join("\t",$zf,"NA","NA",$mouseZFs{$zf})."\n";
}
close ZFS;

open ALLELES, '>', 'PrZFA_allele_codes.mouse.txt';

for my $przfa(sort keys(%pubAlleles)){
  print ALLELES join("\t",$pubAlleles{$przfa},$mouseAlleles{$pubAlleles{$przfa}}->{code})."\n" unless ($seqDone{$mouseAlleles{$pubAlleles{$przfa}}->{seq}}++);
}

for my $przfa(sort keys(%mouseAlleles)){
  print ALLELES join("\t",$przfa,$mouseAlleles{$przfa}->{code})."\n" unless ($seqDone{$mouseAlleles{$przfa}->{seq}}++);
}
close ALLELES;

################################################################################
sub newZF{
  my ($zfFile, $zfRef) = @_;
  open IN, $zfFile;
  while (<IN>){
    chomp;
    my ($name,$seq) = split(/\t/,$_);
    $$zfRef{$name} = $seq;
    $$zfRef{$seq} = $name;
  }
  close IN;
}
