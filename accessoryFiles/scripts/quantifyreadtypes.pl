use strict; 

my ($barcodes_fasta, $prdm9_fa, $write_fastas) = @ARGV;
open IN, $barcodes_fasta;
##'../allBAbarcodes_PacBio.fa';

my %barcode;
my ($BC1,$BC2);

while (<IN>){
	chomp;
	if ($_ =~ /^\>(BC\d+).+(BC\d+)/){
		($BC1,$BC2) = ($1,$2) ;
		next;
	}	
	
	$_ =~ /^.{24}(.{7}).+(.{6}).{3}$/;
	my ($bcSeq1, $bcSeq2) = ($1,$2);
	$barcode{$BC1} = $bcSeq1;
	$barcode{$BC2} = $bcSeq2;
}
close IN;

print join ("\t","BCO","BCC","total_CCS","nNone","nZF","nBCO","nBCC","nZFBCO","nZFBCC","nBC1BC2","nALL","has_multiZF","has_multiBCO","has_multiBCC")."\n";


chomp; 

my $cmd = "cat $prdm9_fa";

$prdm9_fa =~ /(BC\d+).+(BC\d+)/;
my ($BC1,$BC2) = ($1,$2);

my ($zf_lots, $bc1_lots, $bc2_lots) = (0,0,0);
my ($zf, $bc1, $bc2, $bc2noZF, $other) = (0,0,0,0,0);

my ($totCCS,$other,$justZF,$justBC1,$justBC2,$ZFBC1,$ZFBC2,$BC1BC2,$ZFBC1BC2) = split(//,"0"x9);

open ODDSEQ1, '>', join("_",$BC1,$BC2,"noZFnoBarcode.fa");
open ODDSEQ2, '>', join("_",$BC1,$BC2,"noZFwithBarcode.fa");

my $cnt=0;

open my $PIPE, '-|', $cmd;
while (<$PIPE>){
	chomp;
	
	next if ($_ =~ /^\>/);

	my ($count_zfs, $count_bc1, $count_bc2) = (0,0,0);

	my $thisType = 0;

	my $fcBC1 = $barcode{$BC1};
	my $fcBC2 = $barcode{$BC2};

	my $rcBC1 = revcomplement($fcBC1);
	my $rcBC2 = revcomplement($fcBC2);

	$thisType += 1 if ($_ =~ /(TGTGGGC|GGGCTTT|GGGGAGA|GCCCACA|AAAGCCC|TGTGGGC)/);
	$thisType += 2 if ($_ =~ /($fcBC1|$rcBC1)/);
	$thisType += 4 if ($_ =~ /($fcBC2|$rcBC2)/);

	$totCCS++;
	$other++    unless($thisType);
	$justZF++   if ($thisType == 1);
	$justBC1++  if ($thisType == 2);
	$justBC2++  if ($thisType == 4);
	$ZFBC1++    if ($thisType == 3);
	$ZFBC2++    if ($thisType == 5);
	$BC1BC2++   if ($thisType == 6);
	$ZFBC1BC2++ if ($thisType == 7);

	$count_zfs++ while ($_ =~ s/(TGTGGGC|GGGCTTT|GGGGAGA|GCCCACA|AAAGCCC|TGTGGGC)//);
	$count_bc1++ while ($_ =~ s/($barcode{$BC1})//);
	$count_bc2++ while ($_ =~ s/($barcode{$BC2})//);

	$zf_lots++ if ($count_zfs > 5);
	$bc1_lots++ if ($count_bc1 > 5);
	$bc2_lots++ if ($count_bc2 > 5);
	
	print ODDSEQ1 ">".(++$cnt)."\n$_\n" unless ($thisType);
	print ODDSEQ2 ">".(++$cnt)."\n$_\n" if ($thisType =~ /(2|4)/);
}

close ODDSEQ1;
close ODDSEQ2;

print join ("\t",$BC1,$BC2,$totCCS,$other,$justZF,$justBC1,$justBC2,$ZFBC1,$ZFBC2,$BC1BC2,$ZFBC1BC2,$zf_lots,$bc1_lots,$bc2_lots)."\n";

sub revcomplement {
	my $seq = shift;
	
	my $rSeq = reverse $seq;
	
	$rSeq =~ s/G/X/g;
	$rSeq =~ s/A/Y/g;
	$rSeq =~ s/T/Z/g;
	$rSeq =~ s/C/W/g;
	
	$rSeq =~ s/X/C/g;
	$rSeq =~ s/Y/T/g;
	$rSeq =~ s/Z/A/g;
	$rSeq =~ s/W/G/g;
	
	return $rSeq;
}

