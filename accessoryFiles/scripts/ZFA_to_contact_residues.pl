#!/usr/bin/env perl
use strict;

my (%dets, @F);

open IN, $ARGV[0];
while (<IN>){
	chomp; 
	@F=split(/\t/,$_); 
	$dets{$F[0]} = $F[1];
}
close IN;

open I2, $ARGV[1];
while (<I2>){
  chomp;
  my $BS = "";
  my @F = split(/\t/,$_);
  for my $zf( $F[1] =~ m/../g ){
    $BS .= $dets{$zf};
  }
  print $F[0]."\t".$BS."\n";
}

close I2; 
