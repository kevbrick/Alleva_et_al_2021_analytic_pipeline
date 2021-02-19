#!/usr/bin/perl
my $bam = $ARGV[0];
my $prefix = $ARGV[1]?$ARGV[1]."_":"";

my $bamChk = $bam;
$bamChk =~ s/(BC\d+)_rnd\d_([FR])_(BC\d+).+(BC\d+)_rnd\d_([FR])_(BC\d+)//; 

my ($b1,$fr1,$b2,$fr2) = ($1."_".$3,$2,$4."_".$6,$5); 

## Same names and FR/RF orientation
if ($b1 eq $b2 && $fr1 ne $fr2){
	system("mv $bam $prefix$b1.demux.bam");
	system("mv $bam.pbi $prefix$b1.demux.bam.pbi");
	my $xml = $bam; $xml =~ s/\.bam/\.subreadset.xml/;
	system("mv $xml $prefix$b1.subreadset.xml");
}else{
	my $newbam = "bad_barcode_combo_$bam";
	system("mv $bam $newbam");
	system("mv $bam.pbi $newbam.pbi");
	my $xml = $bam;    $xml    =~ s/\.bam/\.subreadset.xml/;
	my $xml = $newbam; $newxml =~ s/\.bam/\.subreadset.xml/;
	system("mv $xml $newxml");
}

