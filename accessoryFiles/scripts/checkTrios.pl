use strict; 

my ($haps,$trioIDs) = @ARGV;

my %haplotypes;
open IN, $haps;
while (<IN>){
  chomp;
  my @F = split(/\t/,$_);
  next if ($_ =~ /\s(NA|Unk|NoZFA)\s/i);
  $haplotypes{$F[0]}->{$F[4]}++;
  $haplotypes{$F[0]}->{$F[5]}++;
  $haplotypes{$F[0]}->{"allele1"} = $F[4];
  $haplotypes{$F[0]}->{"allele2"} = $F[5];
}
close IN;

open TRIO, $trioIDs;
while (<TRIO>){
  chomp;
  my ($kid,$mom,$dad) = split(/\t/,$_);
  next unless ($haplotypes{$kid}->{"allele1"});
  next unless ($haplotypes{$kid}->{"allele2"});
  next unless ($haplotypes{$mom}->{"allele1"});
  next unless ($haplotypes{$mom}->{"allele2"});
  next unless ($haplotypes{$dad}->{"allele1"});
  next unless ($haplotypes{$dad}->{"allele2"});
  my $a1 = $haplotypes{$kid}->{"allele1"};
  my $a2 = $haplotypes{$kid}->{"allele2"};
  
  my ($n1ok,$n2ok,$ok) = (0,0,"NONE");
  
  ($n1ok,$n2ok,$ok) = (1,1,"M1D2") if ($haplotypes{$mom}->{$a1} && $haplotypes{$dad}->{$a2});
  ($n1ok,$n2ok,$ok) = (1,1,"M2D1") if ($haplotypes{$mom}->{$a2} && $haplotypes{$dad}->{$a1});
  
  if ($ok eq "NONE"){
    $ok .= "|M1MOM" if ($haplotypes{$mom}->{$a1});
    $ok .= "|M1DAD" if ($haplotypes{$dad}->{$a1});
    $ok .= "|M2MOM" if ($haplotypes{$mom}->{$a2});
    $ok .= "|M2DAD" if ($haplotypes{$dad}->{$a2});
    
    $n1ok++ if ($haplotypes{$mom}->{$a1} && !($haplotypes{$dad}->{$a2}));
    $n2ok++ if ($haplotypes{$mom}->{$a2} && !($haplotypes{$dad}->{$a1}));
  
  }
  
  print join("\t",$kid,$a1,$a2,
                  $mom,$haplotypes{$mom}->{allele1},$haplotypes{$mom}->{allele2},
                  $mom,$haplotypes{$dad}->{allele1},$haplotypes{$dad}->{allele2},
                  $n1ok,$n2ok,$n1ok+$n2ok,$ok)."\n";
}