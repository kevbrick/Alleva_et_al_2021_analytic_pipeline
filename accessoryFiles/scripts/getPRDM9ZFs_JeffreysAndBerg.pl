#!/usr/bin/perl
use charnames qw[ :full ];

my (%zfs, %zfseq, %zfsrc, %zfseqchk, %allZFalleles, %bloodSpermAlleles);

my $newAlleleName = 0;
my $newAlleleID; 

open ZFS, '>', 'PRDM9_ZFs_Jeffreys_2010.txt';

getZFsFromBerg2010();
getZFsFromJeffreys2012();

getAllelesFromBerg2010();
getAllelesFromJeffreys2012();

getZFsFromBaudat2010();

## OUTPUT ALLELES
open ALLELESBYCODE,     '>', 'humanPRDM9alleles.BergJeffreys.txt';
open ALLELESEQRAW,      '>', 'humanPRDM9alleles.BergJeffreys.withSeq.txt';
open ALLELESFASTA,      '>', 'humanPRDM9alleles.BergJeffreys.fa';
open BLOODSPERMALLELES, '>', 'humanPRDM9alleles.BloodandSpermVariants.txt';

for my $zfs(sort(keys(%allZFalleles))) {
  print ALLELESBYCODE join("\t", $allZFalleles{$zfs}, $zfs)."\n";
  print ALLELESFASTA  '>'.$allZFalleles{$zfs}."===$zfs\n".toSeq($zfs)."\n";
  print ALLELESEQRAW  join("\t",$allZFalleles{$zfs},$zfs,toSeq($zfs,1))."\n";
}

close ALLELESBYCODE;
close ALLELESFASTA;

system('sort -k1,1 humanPRDM9alleles.BergJeffreys.txt -o humanPRDM9alleles.BergJeffreys.txt');

## OUTPUT ZFs
open ALLZFCODES,  '>', 'humanPRDM9ZFcodes.BergJeffreys.txt';
for my $zf(sort(keys(%zfseq))) {
  print ALLZFCODES join("\t",$zf,$zfs{$zf},$zfsrc{$zf},$zfseq{$zf})."\n";
}
close ALLZFCODES;

## OUTPUT BLOOD & SPERM ALLELES
for my $man(sort keys(%bloodSpermAlleles)){
  for my $allele(sort keys(%{$bloodSpermAlleles{$man}})){
    my $alleleName = $allZFalleles{$allele};
    unless ($man =~ /(M\d+[BS]:$alleleName\-|M\d+[BS]:\S+?\-$alleleName)/){
      for my $i (1..$bloodSpermAlleles{$man}->{$allele}){
        print BLOODSPERMALLELES join("\t",$man,$allele,$alleleName)."\n";
      }
    }
  }
}

close BLOODSPERMALLELES;

#########################################################################################################
sub getZFsFromBerg2010{
  system('wget https://static-content.springer.com/esm/art%3A10.1038%2Fng.658/MediaObjects/41588_2010_BFng658_MOESM6_ESM.pdf');

  my $cmd = 'pdftotext -layout -f 13 -l 13 41588_2010_BFng658_MOESM6_ESM.pdf /dev/stdout ';

  open my $PIPE, '-|', $cmd;

  while (<$PIPE>){
    chomp;
    next unless ($_ =~ /^\s+([A-Z])\s+([GATCNX].+)$/);
    my ($name, $seq) = ($1, $2);
    $zfs{"$name"}    = ":$name";
    $zfs{":$name"}   = $name;
    $zfseq{":$name"} = $seq;
    $zfseqchk{$seq} = ":$name";
    $zfsrc{":$name"} = 'BergEtAl2010';
  }

  close ZFs;

}

#########################################################################################################
sub getZFsFromJeffreys2012{
  system('wget https://www.pnas.org/content/pnas/suppl/2012/12/23/1220813110.DCSupplemental/sfig02.pdf');

  my $cmd = 'pdftotext -layout sfig02.pdf /dev/stdout |./replaceUnicode ';

  my $nxt = 'A';

  open my $PIPE, '-|', $cmd;

  while (<$PIPE>){
    chomp;
    next unless ($_ =~ /[GATCNX]{20}/);
    $_ =~ /^(\S+)\s.+?([GATCNX]+)$/;

    my ($init_name, $seq) = ($1,$2);
    my ($name);

    if ($init_name =~ /^[A-Z]$/){
      $name = ":$init_name" ;
    }else{
      if ($init_name =~ /^[a-z|0-9]$/){
        $name = "\|$init_name" ;
      }else{
        $name = "\|".$nxt;
        $nxt = chr(ord($nxt)+1);
      }
    }

    print ZFS    join("\t",$name,$seq)."\n";
    print STDERR join("\t",$init_name,$name,$seq)."\n";

    unless ($zfs{$name}){
      $zfs{"$init_name"} = $name;
      $zfs{"$name"}      = $init_name;
      $zfseq{"$name"}    = $seq;
      $zfseqchk{$seq}   = $name;
      $zfsrc{"$name"}    = 'JeffreysEtAl2012';
    }
  }

  close ZFs;
}

#########################################################################################################
sub getZFsFromBaudat2010{
  system('wget https://science.sciencemag.org/content/sci/suppl/2009/12/30/science.1183439.DC1/Baudat.SOM.pdf');

  system('pdfimages -f 15 -l 16 -png Baudat.SOM.pdf bd');

  system('tesseract bd-000.png bd1');
  system('tesseract bd-001.png bd2');
  system('tesseract bd-002.png bd3');
  system('tesseract bd-003.png bd4');

  system('cat bd1.txt bd2.txt bd3.txt bd4.txt >baudat.txt');

  my $nxt = 'A';
  my (%allele, $name);

  open PIPE, 'baudat.txt';

  while (<PIPE>){
    chomp;
    if ($_ =~ /prdm9-(\S)/i){
      if ($allele{$name}){
        unless ($allZFalleles{$allele{$name}}){
          $allZFalleles{$allele{$name}} = $name;
        }
      }
      $name = $1 ;
      $name = "I" if ($name == 1);
      $name = "baudat_$name";
      print $name."\n";


      next;
    }

    next unless ($_ =~ /TGT.+GGAG$/);

    if (length($_) == 85){
        $_ =~ s/TIT/TTT/;
        $_ =~ s/TI/T/;
    }

    $_ = uc($_);
    $_ =~ s/I/T/g;

    print $_."\n";

    my $seq = $_;

    if ($zfseqchk{$seq}){
      $allele{$name} .= $zfseqchk{$seq};
    }else{
      my $zfname = "b".$nxt;
      $zfseq{$zfname} = $seq;
      $zfs{$zfname} = $zfname;
      $zfseqchk{$seq} = $zfname;
      $zfsrc{$zfname} = 'Baudat_2010';
      $nxt = chr(ord($nxt)+1);
      $allele{$name} .= $zfname;
    }
    my $r = 1;
  }
  if ($allele{$name}){
    unless ($allZFalleles{$allele{$name}}){
      $allZFalleles{$allele{$name}} = $name;
    }
  }
}

#########################################################################################################
sub getAllelesFromBerg2010{
  system('wget https://static-content.springer.com/esm/art%3A10.1038%2Fng.658/MediaObjects/41588_2010_BFng658_MOESM6_ESM.pdf');

  my $cmd = 'pdftotext -layout -f 14 -l 14 41588_2010_BFng658_MOESM6_ESM.pdf /dev/stdout |./replaceUnicode ';

  open my $PIPE, '-|', $cmd;

  while (<$PIPE>){
    chomp;
    next unless ($_ =~ /^\s+([ABCDEL]\d*)\s+\d+\s+(\S+)\s+.+\d+$/);

    my ($name, $allele_seq) = ($1,$2);
    my $allele_code = convertAlleleCode(split(//,$allele_seq));

    $allZFalleles{$allele_code} = $name;

  }
}

#########################################################################################################
sub getAllelesFromJeffreys2012{
  system('wget https://www.pnas.org/content/pnas/suppl/2012/12/23/1220813110.DCSupplemental/sfig04.pdf');

  my $cmd = 'pdftotext -layout sfig04.pdf /dev/stdout |./replaceUnicode ';

  open ALLELES, '>', 'PRDM9_Alleles_Jeffreys_2010.txt';
  open my $PIPE, '-|', $cmd;

  my $newName = 'A';

  while (<$PIPE>){
    chomp;
    if ($_ =~ /Man\s+(\d+),\s+PRDM9\s+(\S+?)\/(\S+),\s(\S+)\smutant/){
      my $type = $4 eq 'sperm'?"S":"B";
      $newAlleleID = "M$1$type:$2-$3";
    }

    next unless ($_ =~ /^\s*(1\-\-|A\S*[BDEFJ])/);

    $_ =~ s/^\s+//;

    my $l = $_;
    if ($_ =~ s/([A-Z|a-z|0-9|\-]{6}\S+)\s+(.+)\s+([A-Z|a-z|0-9|\-]{6}\S+)\s+(.+)//){
      my ($allele1, $dets1, $allele2, $dets2) = ($1,$2,$3,$4);

      storeAllele($allele1, $dets1);
      storeAllele($allele2, $dets2);
      next;
    }

    if ($_ =~ /([A-Z|a-z|0-9|\-]{6}\S+)\s+(.+)/){
      my ($allele1, $dets1) = ($1,$2);
      storeAllele($allele1, $dets1);
      next;
    }

    my $ok = 1;

  }
}

#########################################################################################################
sub storeAllele{
  my ($allele_seq, $allele_dets) = @_;

  $allele_seq =~ s/\-//g;

  my $allele_code = convertAlleleCode(split(//,$allele_seq));

  $bloodSpermAlleles{$newAlleleID}->{$allele_code}++;
  
  return if ($allZFalleles{$allele_code});

  my @F = split(/\s+/,$allele_dets);

  my $name;

  if ($#F == 3){
    $name = $F[3];
  }else{
    if ($#F == 0){
      $name = $F[0];
    }else{
      $F[2] = $F[2] eq '?'?"U":$F[2]; ## For variants of unknown origin
      $F[2] = $F[2] eq "rec"?"R":$F[2]; ## For recombinant variants 
      $name = join(":",$F[2]."v",$F[1],padZeros(++$newAlleleName,4),$newAlleleID);
    }
  }
  $allZFalleles{$allele_code} = $name;
}

#########################################################################################################
sub convertAlleleCode{
  my @aCode = @_;

  my $ret_code;
  for my $ac(@aCode){
    $ret_code .= $zfs{$ac};

    unless ($zfs{$ac}){
      my $r = 1;
    }
  }
  return($ret_code);
}

#########################################################################################################
sub toSeq{
  my ($zfcode, $noEL) = @_;

  my $retSeq;

  while ($zfcode =~ s/^(..)//s){
    $retSeq .= $zfseq{$1}.($noEL?"":"\n");
  }

  chomp $retSeq unless ($noEL);

  return $retSeq;

}

#########################################################################################################
sub padZeros{
  my ($num, $len) = @_;

  $num = "0$num" while (length($num) < $len);
  return $num;
}
