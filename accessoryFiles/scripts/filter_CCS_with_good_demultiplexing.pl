#!/usr/bin/perl
use strict;
use Getopt::Long;

GetOptions ('p=s'      => \(my $popFile),
            'n=s'      => \(my $nameFile),
						'e=s'      => \(my $exptName));

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

## GET NAME TO BARCODE TABLE#####################################
my %names;
my %edate;

open NM, $nameFile;

while (<NM>){
  chomp;
	next unless ($_ =~ /($exptName)/);
  my ($name,$bc2,$bc1,$expt,$edate)  = split(/,/,$_);

  $names{$bc1."_".$bc2} = $name;
  $edate{$expt} = $edate;
}

close NM;

#################################################################
## Loop through all demultiplexed BAMs
opendir(DIR,'.');

while (my $f = readdir(DIR)){

	next unless ($f =~ /BC.+BC.+\.bam$/);

	my $bam    = $f;
	my $bamChk = $bam;

	$bamChk =~ s/BC(\d+)_([FR])_BC(\d+).+BC(\d+)_([FR])_BC(\d+)//;

	my ($bc1A,$fr1,$bc1B,$bc2A,$fr2,$bc2B) = ($1,$2,$3,$4,$5,$6);

	my ($b1,$b2) = ($bc1A."_".$bc1B, $bc2A."_".$bc2B);

	## DMX good IF: Same names and FR/RF orientation
	if ($b1 eq $b2 && $fr1 ne $fr2){

		my $name = $names{$b1};

		if (!$pop{$name}){
      $pop{$name}->{popid} = 'other';
      $pop{$name}->{pop}   = 'OTH';
      $pop{$name}->{sex}   = 'MALE';
    }

    $name     = join("_",$name,$pop{$name}->{pop},$bc1A,$bc1B,uc($pop{$name}->{sex}) eq "MALE"?"M":"F");

		my $xml = $bam;
		$xml =~ s/\.bam/\.subreadset.xml/;

		my $newbam    = $name."_".$exptName.".demux.bam";
		my $newbampbi = $name."_".$exptName.".demux.bam.pbi";
		my $newxml    = $name."_".$exptName.".subreadset.xml";
    my $newfa     = $name."_".$exptName.".raw.fa";

		system("mv $bam $newbam");
		system("mv $bam.pbi $newbampbi");
   	system("mv $xml ".$name."_".$exptName.".subreadset.xml");

		open FA, '>', $newfa;

		my $seqCnt;
		open my $IN, '-|', "samtools view $newbam |cut -f10";
		while (<$IN>){
			chomp;
			print FA ">".$name."_".$exptName."|".(++$seqCnt)."\n";
			print FA to84($_);
		}
		close FA; close $IN;
	}else{
		my $newbam = $bam; $newbam =~ s/demux.bam/demuxBAD.bam/;
		system("mv $bam $newbam");
		system("mv $bam.pbi $newbam.pbi");
		my $xml = $bam;    $xml       =~ s/\.bam/\.subreadset.xml/;
		my $xml = $newbam; my $newxml =~ s/\.bam/\.subreadset.xml/;
		system("mv $xml $newxml");
	}
}

sub to84{
  my $inLine = shift;
  my $lnRet;
  while ($inLine =~ s/^(.{84})//){
    $lnRet .= $1."\n";
  }
  $lnRet .= $1."\n";
  return($lnRet);
}
