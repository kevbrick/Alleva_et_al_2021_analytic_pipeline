#!/usr/bin/perl
use strict;

my $folder = $ARGV[0];

my (%gt, %newZFs, %newAlleles, %newZFDets, %newAlleleDets);
my ($newZFCounter, $newAlleleCounter) = (65,1);

getHaplotypes($folder,\%gt);

print join("\t",'id',
                'pop',
	            	'diploid',
                'homhet',
                'allele_1',
                'allele_2',
                'allele_1_code',
                'allele_2_code',
                'size_1',
                'size_2',
                'seq_1',
                'seq_2')."\n";

for my $id (sort keys(%gt)){
		print join("\t",$id,
	    $gt{$id}->{'pop'},
		  $gt{$id}->{'diploid'},
		  $gt{$id}->{'homhet'},
		  $gt{$id}->{'allele_1'},
		  $gt{$id}->{'allele_2'},
		  $gt{$id}->{'allele_1_code'},
		  $gt{$id}->{'allele_2_code'},
		  $gt{$id}->{'size_1'},
		  $gt{$id}->{'size_2'},
		  $gt{$id}->{'seq_1'},
		  $gt{$id}->{'seq_2'})."\n";
}

open NEWZFS, '>', 'newPRDM9ZFs.txt';
open NEWALLELES, '>', 'newPRDM9Alleles.txt';

for my $z(keys %newZFDets){
	print NEWZFS join("\t",$z,$newZFDets{$z}->{cnt},$newZFDets{$z}->{seq})."\n";
}

for my $z(keys %newAlleleDets){
	print NEWALLELES join("\t",$z,$newAlleleDets{$z}->{cnt},$newAlleleDets{$z}->{code})."\n";
}

close NEWZFS;
close NEWALLELES;

sub orNA{
  my $v = shift;

  return ($v?$v:"NA");

}

sub printAll{
	my ($txt,@arr) = @_;

	my $ret;
	for my $x(@arr){$ret .= $txt."_$x"."\t"};
	chop $ret;
	return($ret);
}

sub printAllVal{
	my ($txt,$ref,@arr) = @_;

	my $ret;
	for my $x(@arr){$ret .= orNA($$ref->{$txt."_$x"})."\t"};
	chop $ret;
	return($ret);
}

#############################################################
sub getHaplotypes{
  my ($hapFolder,$genotypes) = @_;

  open my $PIPE, '-|', 'find '.$hapFolder.' -name "*haplotypes.txt" |xargs cat';

  my $c;

  while (<$PIPE>){
  	chomp $_;
  	my ($uid,$allele,$alleleCode,$size,$count,$seq) = split(/\t/,$_);
  	$uid =~ /^(\S+)_(\S+)/;
		my ($id, $pop) = ($1,$2);

		#next unless ($id && $pop && $allele && $alleleCode && $size && $count && $seq);
    next unless ($id && $pop && $allele);

  	if ($allele eq 'new'){
	  	if ($seq){
		  	($allele, $alleleCode) = dealwithNewZFsandAlleles($alleleCode,$seq) ;
		  }else{
  			($allele, $alleleCode) = ('NA','NANANANANANA');
				$size  = $size?$size:"NA";
				$count = $count?$count:"NA";
				$seq   = $seq?$seq:"NA";
		  }
	  }

		$$genotypes{$id}->{'pop'} = $pop;

		if ($$genotypes{$id}->{'allele_1'}){
			## assure alleles remain in alphabetical order
			my @aOrd = sort($$genotypes{$id}->{'allele_1'},$allele);
			if ($$genotypes{$id}->{'allele_1'} eq $aOrd[0]){
			  $$genotypes{$id}->{'allele_2'}      = $allele;
			  $$genotypes{$id}->{'allele_2_code'} = $alleleCode;
			  $$genotypes{$id}->{'size_2'}        = $size;
	  	  $$genotypes{$id}->{'seq_2'}         = $seq;
			}else{
			  $$genotypes{$id}->{'allele_2'}      = $$genotypes{$id}->{'allele_1'};
			  $$genotypes{$id}->{'allele_2_code'} = $$genotypes{$id}->{'allele_1_code'};
			  $$genotypes{$id}->{'size_2'}        = $$genotypes{$id}->{'size_1'};
			  $$genotypes{$id}->{'seq_2'}         = $$genotypes{$id}->{'seq_1'};

			  $$genotypes{$id}->{'allele_1'}      = $allele;
			  $$genotypes{$id}->{'allele_1_code'} = $alleleCode;
			  $$genotypes{$id}->{'size_1'}        = $size;
	  		$$genotypes{$id}->{'seq_1'}         = $seq;
			}

			$$genotypes{$id}->{'diploid'}       = $$genotypes{$id}->{'allele_1'}."/".$$genotypes{$id}->{'allele_2'};
	  	$$genotypes{$id}->{'homhet'}        = 'het';
	  }else{
		  $$genotypes{$id}->{'allele_1'}      = $allele;
		  $$genotypes{$id}->{'allele_1_code'} = $alleleCode;
		  $$genotypes{$id}->{'size_1'}        = $size;
	  	$$genotypes{$id}->{'seq_1'}         = $seq;
		}
  }

  for my $id (keys %{$genotypes}){
	  unless ($$genotypes{$id}->{'allele_2'}){
  		$$genotypes{$id}->{'allele_2'}      = $$genotypes{$id}->{'allele_1'};
		  $$genotypes{$id}->{'allele_2_code'} = $$genotypes{$id}->{'allele_1_code'};
		  $$genotypes{$id}->{'size_2'}        = $$genotypes{$id}->{'size_1'};
		  $$genotypes{$id}->{'seq_2'}         = $$genotypes{$id}->{'seq_1'};

		  $$genotypes{$id}->{'diploid'}       = $$genotypes{$id}->{'allele_1'}."/".$$genotypes{$id}->{'allele_1'};
		  $$genotypes{$id}->{'homhet'}        = 'hom';
	  }
  }
  close $PIPE;
}

#################################################################
sub newZF{
	my $s = shift;
	my @ret;
	while ($s =~ s/^(.{84})//){push @ret, $1};
	return (@ret);
}

#################################################################
sub dealwithNewZFsandAlleles{
	my ($aCode, $zfseq) = @_;
	my (@arrA,@arrZ);
	while ($aCode =~ s/^(.{2})//){push @arrA, $1};
	while ($zfseq =~ s/^(.{84})//){push @arrZ, $1};

	for my $i(0..$#arrA){
		if ($arrA[$i] eq 'NU'){
			if ($newZFs{$arrZ[$i]}){
				$arrA[$i] = $newZFs{$arrZ[$i]};
				$newZFDets{$arrA[$i]}->{count}++;
			}else{
				$arrA[$i] = '!'.chr($newZFCounter++);
				$newZFs{$arrZ[$i]} = $arrA[$i];
				$newZFDets{$arrA[$i]}->{count}++;
				$newZFDets{$arrA[$i]}->{seq} = $arrZ[$i];
			}
		}
	}

	my $retCode = join("",@arrA);
	my $retAllele;
	if ($newAlleles{$retCode}){
		$retAllele = $newAlleles{$retCode};
		$newAlleleDets{$retAllele}->{cnt}++;
	}else{
		$retAllele = "M".($newAlleleCounter++);
		$newAlleles{$retCode} = $retAllele;
		$newAlleleDets{$retAllele}->{cnt}++;
		$newAlleleDets{$retAllele}->{code} = $retCode;
	}

	return ($retAllele,$retCode);
}
