#!/usr/bin/perl
## Script to combine the barcode1 and barcode2 sequences from Ben's
## amplicaon sequencing strategy
use strict;

my (%fa1,$name);

open FA1, 'barcode_96.fasta';

while (<FA1>){
  chomp;
  if ($_ =~ /^\>(\S+)/){
    $name = $1; next;
  }else{
    ## Add FWD and REV sequences that are present in ONT kit
    ## Each barcode will have both
    $fa1{$name.":R"} = $_."TTAACCTACTTGCCTGTCGCTCTATCTTCGGCGTCTACTTGGGTGTTTAACCT";
    $fa1{$name.":F"} = $_."TTAACCTTTCTGTTGGTGCTGATATTGCGGCGTCTGCTTGGGTGTTTAACCT";
    $name = '';
  }
}

close FA1;

## Here, we make all possible bc1-bc2 combinations
$name = '';
open FA2, 'barcodes_custom.12bp.fa';
while (<FA2>){
  chomp;
  if ($_ =~ /^\>(\S+)/){
    $name = $1; next;
  }else{
    for my $faName1(sort keys(%fa1)){
      my $seq = $fa1{$faName1}.$_;
      print ">".$faName1.":".$name."\n";
      print "$seq\n";
    }
  }
}
