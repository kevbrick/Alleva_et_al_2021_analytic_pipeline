
use strict; 
use Math::Round; 

open ZF, '/data/RDCO/code/isoseq3/dataOut/parseZFs2/humanPRDM9ZFcodes.BergJeffreys.txt';

my %zfs; 
my %zfCounts; 

while (<ZF>){
	chomp;
	my @F = split(/\t/,$_);
	$zfs{$F[0]} = $F[3];
	$zfs{$F[3]} = $F[0];
}

close ZF;

open PRDM9, '/data/RDCO/code/isoseq3/dataOut/parseZFs2/humanPRDM9alleles.BergJeffreys.txt';

my %prdm9; 
my %prdm9Counts; 

while (<PRDM9>){
	chomp;
	my @F = split(/\t/,$_);
#	$F[1] =~ s/\:J\s*$/_0/;
	$prdm9{$F[1]} = $F[0];

}

close PRDM9;

my %ontZF;
open ONT, '/data/RDCO/mini/PRDM9/100319_prdm9_758samples/updated_prdm9_alleles_1d2.txt';
while(<ONT>){
	chomp;
	my ($name, $zfarray_code, $zfarray_len, $zfarray_seq, $zfarray_n) = split(/\t/,$_);

	$name =~ s/BC99/BC105/;
	$name =~ s/BC100/BC106/;

	if ($ontZF{$name}){
		$ontZF{$name."a"} = $zfarray_seq;
	}else{
		$ontZF{$name} = $zfarray_seq;
	}
}

my %pops;
$pops{"BC".$_} = 'Yoruba(Nigeria)' foreach(1..15);
$pops{"BC".$_} = 'Lhuya(Kenya)'    foreach(16..30);
$pops{"BC".$_} = 'Italian'         foreach(31..45);
$pops{"BC".$_} = 'Finnish'         foreach(46..50);
$pops{"BC".$_} = 'Peruvian'        foreach(59..67);
$pops{"BC".$_} = 'Pakistani'       foreach(68..81);
$pops{"BC".$_} = 'Han(Chinese)'    foreach(82..96);

my $nuChr = chr(48);
opendir(DIR,'.');

my %okCnt;
while (my $chk_file = readdir(DIR)){

	next unless ($chk_file =~ /\.raw$/);

	$chk_file =~ /(BC\d+).+?(BC\d+)/;

	my ($bc1,$bc2) = ($1,$2);

	my ($alleleCount, $consensus, $seq, $a1_seq, $a2_seq, $a1_zfcode, $a2_zfcode);

	open RAW, $chk_file; 
	while (<RAW>){
		chomp;
		my @F = split(/\t/,$_);
		
		if ($F[0] == 1 && ++$alleleCount == 2){
			$a1_zfcode = $consensus;
			$a1_seq = $seq;
			$consensus = '';
			$seq = '';
		}

		if ($zfs{$F[1]}){
			$consensus .= $zfs{$F[1]};
			$seq       .= $F[1];
			$zfCounts{$zfs{$F[1]}}++;
		}else{	
			my $newZF    = '_'.$nuChr;			
			$zfs{$F[1]}  = $newZF;
			$zfs{$newZF} = $F[1];
			$zfCounts{$newZF}++;

			$nuChr       = chr(ord($nuChr)+1);
			$nuChr       = chr(65) if (ord($nuChr) == 57);
			$nuChr       = chr(97) if (ord($nuChr) == 91);
			$consensus  .= $zfs{$F[1]};
                        $seq        .= $F[1];
		}
		
	}

	if ($a1_zfcode){
		$a2_zfcode = $consensus;
		$a2_seq    = $seq;
	}else{
		$a1_zfcode = $consensus;
		$a2_zfcode = $consensus;
		$a1_seq    = $seq;
		$a2_seq    = $seq;
	}

	my ($a1_name, $a2_name, $a1_len, $a2_len);
	my $homhet = 'unknown';


	if ($a1_seq){
		$a1_len  = length($a1_zfcode)/2;
		$a1_name = $prdm9{$a1_zfcode}?$prdm9{$a1_zfcode}:"new";
		$homhet  = 'hom';	

		if ($a2_zfcode){
			$a2_len  = length($a2_zfcode)/2;
			$a2_name = $prdm9{$a2_zfcode}?$prdm9{$a2_zfcode}:"new";
			$homhet  = 'het' if ($a1_zfcode ne $a2_zfcode);
		}	
	}else{
		$a1_zfcode = 'NA';
		$a1_seq    = 'NA';
		$a1_len    = 'NA';
		$a1_name   = 'NA';
		$a2_seq    = 'NA';
		$a2_zfcode = 'NA';
		$a2_len    = 'NA';
		$a2_name   = 'NA';
		$okCnt{NA}++;
		
	}

	my ($allele1InONT, $allele2InONT) =("N","N");

	my $ontName = $bc1."_".$bc2;

	open ONTCOMP, '>', $ontName.".PB_v_ONT.ZFarrays.fa";

	print ONTCOMP ">PB_allele1\n$a1_seq\n";
	print ONTCOMP ">PB_allele2\n$a2_seq\n" if ($a1_seq ne $a2_seq);
	print ONTCOMP ">ONT_allele1\n$ontZF{$ontName}\n";
	print ONTCOMP ">ONT_allele2\n".$ontZF{$ontName."a"}."\n" if ($ontZF{$ontName."a"});
	close ONTCOMP;

 	$allele1InONT = "Y" if ($ontZF{$ontName}                             eq $a1_seq);
 	$allele1InONT = "Y" if ($ontZF{$ontName."a"} && $ontZF{$ontName."a"} eq $a1_seq);
 	$allele2InONT = "Y" if ($ontZF{$ontName}                             eq $a2_seq);
 	$allele2InONT = "Y" if ($ontZF{$ontName."a"} && $ontZF{$ontName."a"} eq $a2_seq);

	$allele1InONT = "A" unless ($ontZF{$ontName});
	$allele2InONT = "A" unless ($ontZF{$ontName});

	my $okType = $allele1InONT.$allele2InONT;

	my $type;
	$type = 'hom'         if ($homhet eq 'hom');
	$type = 'tough_het'   if ($homhet eq 'het' && $a1_len == $a2_len); 
	$type = 'simple_het'  if ($homhet eq 'het' && $a1_len != $a2_len);

	my $homhetONT = ($ontZF{$ontName} && $ontZF{$ontName."a"})?"ONThet":"ONThom";

	if ($okType eq "YY"){
		$okCnt{$type}->{OK}++;
	}else{
		$okCnt{$type}->{BAD}++;
	}
 	
	print join("\t", $bc1, $bc2, $a1_zfcode, $a2_zfcode, $a1_name, $a2_name, $a1_len, $a2_len, $homhet, $type, $allele1InONT, $allele2InONT, $homhetONT, join("/",sort($a1_name,$a2_name)), $pops{$bc1})."\n";

}

print join("\t","type","ok","bad","okpc","badpc","tot")."\n";
for my $t('hom','simple_het','tough_het'){
	$okCnt{$t}->{OK}  = $okCnt{$t}->{OK}?$okCnt{$t}->{OK}:0;
	$okCnt{$t}->{BAD}  = $okCnt{$t}->{BAD}?$okCnt{$t}->{BAD}:0;

	$okCnt{$t}->{TOT} = $okCnt{$t}->{OK} + $okCnt{$t}->{BAD};

	print join("\t", $t,$okCnt{$t}->{OK},$okCnt{$t}->{BAD},round($okCnt{$t}->{OK}/$okCnt{$t}->{TOT}*100),round($okCnt{$t}->{BAD}/$okCnt{$t}->{TOT}*100),$okCnt{$t}->{TOT})."\n";
}

print join("\t", "unknown",$okCnt{NA},0,100,0,$okCnt{NA})."\n";

