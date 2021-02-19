#!/usr/bin/perl
use strict;

open IN, $ARGV[0];
my $name = $ARGV[1]?$ARGV[1]."_":"";

my ($motif,$motifLen);
my ($A,$C,$G,$T);

while (<IN>){
  if ($_ =~ /Motif\s+(\S+)\s.+position-specific\s+prob/){$motif = $1; next}
  if ($_ =~ /letter-probability.+w\=\s+(\d+)/){$motifLen = $1; next}

  next if ($_ =~ /^[\s\-]+$/);

  if ($motif){
    my $r = 1;
  }

  if ($motifLen){
    if ($_ =~ /^\s*[01]\.\d+/){
      chomp;
      $_ =~ s/^\s+//;
      $_ =~ s/\s+(\S)/\t$1/g;
      my @F = split(/\t/,$_);
      $motifLen--;
      $A .= $F[0]."\t";
      $C .= $F[1]."\t";
      $G .= $F[2]."\t";
      $T .= $F[3]."\t";
    }
  }else{
    if ($motif){
      open PFM, '>', $name."MEME_motif_$motif.pfm";
      print PFM "A\t$A\n";
      print PFM "C\t$C\n";
      print PFM "G\t$G\n";
      print PFM "T\t$T\n";
      close PFM;
      ($motif,$motifLen) = ('',0);
      ($A,$C,$G,$T) = ();
    }
  }
}
