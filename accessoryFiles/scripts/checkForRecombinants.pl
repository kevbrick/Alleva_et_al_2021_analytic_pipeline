use strict;
use Getopt::Long;

GetOptions ('allele=s' => \(my $alleles),
            'all=s'    => \(my $all_alleles_file = 'PrZFA_alleles.details.txt'));

my @alleles_to_check;

## Allow input as a single allele, a comma separated list, or a file
if (-e $alleles){
  open ALL, $alleles;
  chomp(@alleles_to_check = <ALL>); 
  close ALL;
}else{
  if ($alleles =~ /,/){
    @alleles_to_check = split(/,/,$alleles);
  }else{
    push @alleles_to_check, $alleles; 
  }
}

my @alleles;
open IN, $all_alleles_file;
while (<IN>){
  my @F = split(/\t/,$_);
  push @alleles, $F[1] unless ($_ =~  /^(#|ID)/);
}

print join("\t","ID","P1ID","P2ID","type","switch_string","source_string","num_switches")."\n";

for my $allele_name(@alleles_to_check){
  for my $i(0..$#alleles){
    my $p1 = $alleles[$i];
    unless ($p1 eq $allele_name){
      
      my $out = `perl checkPrdm9BiParentalRecombinants.pl --n $p1,$p1,$allele_name --mode 1 --all $all_alleles_file`;
      chomp $out;
      print $out."\n";    
      
      for my $j(($i+1)..$#alleles){
        my $p2 = $alleles[$j];
        my $aStr = "$p1,$p2,$allele_name";
        unless ($p1 eq $allele_name || $p2 eq $allele_name){
          my $out = `perl checkPrdm9BiParentalRecombinants.pl --n $aStr --mode 2 --all $all_alleles_file`;
          chomp $out;
          print $out."\n";    
        }
      } 
    }
  }
}
