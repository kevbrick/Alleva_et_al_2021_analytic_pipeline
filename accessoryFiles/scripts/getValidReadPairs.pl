#!/usr/bin/perl
use strict;
use Getopt::Long;

GetOptions ('p=s'      => \(my $popFile),
            'n=s'      => \(my $nameFile),
            'bc=s'     => \(my $barcodeFile),
            'pairs=s'  => \(my $pairFile),
            'out=s'    => \(my $outFile));

#my ($popFile, $nameFile, $barcodeFile, $pairFile, $outFile) = @ARGV;

## GET INDIVIDUAL DETAILS########################################
my %pop;

open POP, $popFile;

while (<POP>){
  chomp;
  my ($pop,$popid,$id,$sex)  = split(/\t/,$_);

  $pop{$id}->{popid} = $popid;
  $pop{$id}->{pop}   = $pop;
  $pop{$id}->{sex}   = $sex;
}

close POP;

## GET NAMES #####################################################
my %names;
my %edate;

open NM, $nameFile;

while (<NM>){
  chomp;
  my ($name,$bc2,$bc1,$expt,$edate)  = split(/,/,$_);

  $names{$bc1."_".$bc2} = $name;
  $edate{$expt} = $edate;
}

close NM;

## GET BARCODES ##################################################
my %ids;

open BC, $barcodeFile;

while (<BC>){
  chomp;
  my ($id,$bc1,$bc2)     = split(/\t/,$_);
  #$barcodes{"$bc1_$bc2"} = $id;
  $ids{$id}              = $bc1."_".$bc2;
}

close BC;

##################################################################
open IDOUT  , '>', $outFile;

##################################################################
open PAIR, $pairFile;

my $lineCount;

while (<PAIR>){
  chomp;

  my ($id1,$id2)     = split(/\t/,$_);

  unless ($lineCount++){
      print IDOUT join("\t",'id1','id2','bc11','bc12','bc21','bc22','name','individual','pop','popid','sex')."\n";
  }

  ## Get barcodes
  my ($bc1A,$bc1B,$bc2A,$bc2B);
  if ($ids{$id1}){
    ($bc1A,$bc1B) = split(/_/,$ids{$id1});
  }else{
    ($bc1A,$bc1B) = ("unk","unk");
  }

  if ($ids{$id2}){
    ($bc2A,$bc2B) = split(/_/,$ids{$id2});
  }else{
    ($bc2A,$bc2B) = ("unk","unk");
  }

  ##
  my $name = "unk";

  if ($ids{$id1} && ($ids{$id1} eq $ids{$id2})){
    $name = $names{$ids{$id1}} ;
  };

  my $idString;

  if ($name && $name ne 'unk'){
    if (!$pop{$name}){
      $pop{$name}->{popid} = 'other';
      $pop{$name}->{pop}   = 'OTH';
      $pop{$name}->{sex}   = 'MALE';
    }
    $idString = join("\t",$name,$pop{$name}->{pop},$pop{$name}->{popid},$pop{$name}->{sex});
    $name     = join("_",$name,$pop{$name}->{pop},$bc1A,$bc1B,uc($pop{$name}->{sex}) eq "MALE"?"M":"F");
  }else{
    $idString = join("\t","UNK","UNK","UNK","UNK");
    $name = "UNK";
  }

  $name =~ s/\//\-/;

  print IDOUT   join("\t",$id1,$id2,$bc1A,$bc1B,$bc2A,$bc2B,$name,$idString)."\n";

}

close PAIR; close IDOUT;
