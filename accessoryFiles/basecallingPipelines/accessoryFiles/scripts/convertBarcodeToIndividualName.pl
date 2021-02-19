use strict; 

#!/usr/bin/perl
use strict;
use Getopt::Long;

GetOptions ('p=s'      => \(my $popFile),
            'n=s'      => \(my $nameFile),
            'out=s'    => \(my $outFile));

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

#################################################################
#my $expt = $ARGV[0];

my %nm;

open IN, $nameFile;

while (<IN>){
	chomp; 
	my ($id,$bc2,$bc1,$sample,$date) = split(/,/,$_);
	next unless ($sample eq $expt);

	$nm{"BC$bc1"."_BC$bc2"} = $id; 
}
close IN; 

opendir(DIR,'.');
while (my $f = readdir(DIR)){
	chomp $f;
	if ($f =~ /^(.*)(BC\d+_BC\d+).+(.demux.bam|subreadset.xml)/){
		system("mv $f ".join("",$nm{$2},".expt_",$expt,"$3")) if ($nm{$2});	
	}

	if ($f =~ /^(.*)_((NA|HG)\S+?).exp.+(.demux.bam|subreadset.xml)/){
		system("mv $f ".join("",$2,".expt_",$expt,"$3")) if ($nm{$2});	
	}
}
