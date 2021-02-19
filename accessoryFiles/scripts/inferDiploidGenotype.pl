use strict;

open ZF, 'humanPRDM9ZFcodes.BergJeffreys.txt';

my %zfs;
my %zfCounts;

while (<ZF>){
	chomp;
	my @F = split(/\t/,$_);
	$zfs{$F[0]} = $F[3];
	$zfs{$F[3]} = $F[0];
}

close ZF;

open PRDM9, 'humanPRDM9alleles.BergJeffreys.txt';

my %prdm9;
my %prdm9Counts;

while (<PRDM9>){
	chomp;
	my @F = split(/\t/,$_);
#	$F[1] =~ s/\:J\s*$/_0/;
	$prdm9{$F[1]} = $F[0];

}

close PRDM9;

my $nuChr = chr(48);
opendir(DIR,'.');

while (my $chk_file = readdir(DIR)){

	next unless ($chk_file =~ /.raw$/);

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
	}


	print join("\t", $bc1, $bc2, $a1_zfcode, $a2_zfcode, $a1_name, $a2_name, $a1_len, $a2_len, $homhet, $a1_seq, $a2_seq)."\n";

}
