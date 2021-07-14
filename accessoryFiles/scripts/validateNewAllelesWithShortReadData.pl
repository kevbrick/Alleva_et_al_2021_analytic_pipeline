#!/usr/bin/perl
use strict; 
use Getopt::Long;

my $skipCheck = 1;  
 
GetOptions ('check'  => sub { $skipCheck = 0 },
            'dets=s' => \(my $dets_file = 'PrZFA_ZFs.details.txt'),
            'haps=s' => \(my $haplotypes_file = 'prdm9_haplotypes.final.tab'),
            'noprefix+' => \(my $noPrefix),
            'ftpprefix=s' => \(my $ftpPrefix = 'ftp:/ftp.1000genomes.ebi.ac.uk/vol1/'),
            'tree=s' => \(my $ftpTree = 'ftp_1KG.current.tree.June232021.txt'));

################################################################################
my %URLs; get1KGURLs();

################################################################################
my $zfFasta = 'PrZFA_ZFs.fasta';

open ZFS, $dets_file;
open ZFFA, '>', $zfFasta;

while (<ZFS>){
  chomp; 
  next unless ($_ =~ /^\S\S\s/);
  
  my @F = split(/\t/,$_);
  print ZFFA '>'.$F[0]."\n".$F[3]."\n";
}
close ZFS;
close ZFFA;
################################################################################
open ALL, $haplotypes_file ;

my $lineNum; 

while (<ALL>){
  chomp; 
  next if ($_ =~ /homhet\s+allele_1\s+allele_2/);
  
  my ($id,$pop,$diploid,$homhet,$allele_1,$allele_2,$allele_1_code,$allele_2_code,$size_1,$size_2,$seq_1,$seq_2,$nseqs,$nzfseqs,$is_child) = split(/\t/,$_);
  
  next unless ($allele_1 =~ /^M\d+$/ || $allele_2 =~ /^M\d+$/);
  
  my %zftypes;
  uniqueZFs(\%zftypes,$allele_1_code,$allele_2_code,$seq_1,$seq_2,$allele_1,$allele_2);
  #my @onlyA2 = uniqueZFs(\$zftypes,$allele_2_code,$allele_1_code,$seq_2,$seq_1);
  
  unless ($URLs{$id}){
    print STDERR "No data for $id !!\n" ;
    next;
  }else{
    print STDERR $id.": ".($URLs{$id}->{exome}?$URLs{$id}->{exome}:"NOEXOME")."\n";
  }
}
close ALL;

###################################################
open HAPS, $haplotypes_file ;

my $lineNum; 

while (<HAPS>){
  chomp; 
  next if ($_ =~ /homhet\s+allele_1\s+allele_2/);
  
  my ($id,$pop,$diploid,$homhet,$allele_1,$allele_2,$allele_1_code,$allele_2_code,$size_1,$size_2,$seq_1,$seq_2,$nseqs,$nzfseqs,$is_child) = split(/\t/,$_);
  
  next unless ($allele_1 =~ /^M\d+$/ || $allele_2 =~ /^M\d+$/);
  
  my %zftypes;
  uniqueZFs(\%zftypes,$allele_1_code,$allele_2_code,$seq_1,$seq_2,$allele_1,$allele_2);

  unless ($URLs{$id}){
    print STDERR "No data for $id !!\n" ;
    next;
  }
  
  print STDERR "$id : ".($URLs{$id}->{exome}?'EXOME':($URLs{$id}->{high_coverage}?'HI-COV':'LO-COV'))."\n";
  my $url = $URLs{$id}->{exome}?$URLs{$id}->{exome}:($URLs{$id}->{high_coverage}?$URLs{$id}->{high_coverage}:$URLs{$id}->{low_coverage});  
  
  my $header   = $id.'.1KG.header';
  my $bam1KG   = $id.'.1KG.bam';
  my $cram1KG  = $id.'.1KG.cram';
  my $bamAln   = $id.'.aln.bam';
  my $bamAlnS  = $id.'.aln.sorted.bam';
  my $samFinal = $id.'.final_filtered.sam';
  my $bamFinal = $id.'.final_filtered.bam';
  my $fq1KG    = $id.'.1KG.fq';
  
  if ($skipCheck || !(-e $bam1KG)){
    my ($chrom,$prdm9Pos);
    system("samtools view -H $url >$header");
    $chrom    = "chr5"                       if (`grep SN:chr5 $header`);
    $chrom    = "5"                          if (`grep SN:5    $header`);
    # $prdm9Pos = $chrom.":23525000-235320000" if (`grep GRCh37  $header`);
    # $prdm9Pos = $chrom.":23525000-235320000" if (`grep NCBI37  $header`);
    # $prdm9Pos = $chrom.":23525000-235320000" if (`grep GRCh38  $header`);
    $prdm9Pos = $chrom.":23525000-23532000";
    
    system('samtools view -bh '.$url.' '.$prdm9Pos.' >'.$bam1KG) ;
  }
  system('bedtools bamtofastq -i '.$bam1KG.' -fq '.$fq1KG)                            if ($skipCheck || !(-e $fq1KG));
  system('minimap2 -x sr -a '.$zfFasta.' '.$fq1KG.' |samtools view -Shb - >'.$bamAln) if ($skipCheck || !(-e $bamAln));
  system('samtools sort '.$bamAln.' >'.$bamAlnS)                                      if ($skipCheck || !(-e $bamAlnS));
  system('samtools index '.$bamAlnS);
  
  open SAM, '>', $samFinal;
  open my $PIPE, '-|', 'samtools view -h -F 4 '.$bamAlnS;
  while (<$PIPE>){
    chomp;
    my @F = split(/\t/,$_);
    if ($_ =~ /^@/ || ($F[5] !~ /GATC/ && $F[5] =~ /(7[3-9]|8[01234])M/)){
      print SAM $_."\n";
    }
  }
  close $PIPE;
  close SAM;
  
  system('samtools view -Shb '.$samFinal.' >'.$bamFinal);
  system('samtools index '.$bamFinal);
  
  open my $PIPE2, '-|', 'samtools idxstats '.$bamFinal;
  while (<$PIPE2>){
    chomp;
    my ($zfID, $zfLen, $coverage) = split(/\t/,$_);
    next if ($zfID eq "*");
    print join("\t",$id,
                    $diploid,
                    $zfID,
                    ($zftypes{$zfID}->{type}?$zftypes{$zfID}->{type}:"None"),
                    $coverage,
                    ($zftypes{$zfID}->{count}?$zftypes{$zfID}->{count}:0),
                    ($zftypes{$zfID}->{count}?$coverage/$zftypes{$zfID}->{count}:$coverage))."\n";
  }
  close $PIPE2;
}
close HAPS;

################################################################################
sub uniqueZFs{
  my ($zftypes,$codeA,$codeB,$seqA,$seqB,$alleleA,$alleleB) = @_;

  my %zfsA;
  my %zfsB;

  ## Get ZFs
  my @arrayA = ( $codeA =~ m/../g );
  for my $zf(@arrayA) {$zfsA{$zf}++; $$zftypes{$zf}->{count}++}
    
  my @arrayB = ( $codeB =~ m/../g );
  for my $zf(@arrayB) {$zfsB{$zf}++; $$zftypes{$zf}->{count}++}
  
  my @retArray;
  for my $zf(@arrayA) {$$zftypes{$zf}->{type} = $zfsB{$zf}?"Both":$alleleA}
  for my $zf(@arrayB) {$$zftypes{$zf}->{type} = $zfsA{$zf}?"Both":$alleleB}
}

################################################################################
sub get1KGURLs{
  
  if (!(-e $ftpTree)){
    system('wget -O '.$ftpTree.' ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/current.tree')
  }
  
  open TREE1KG, $ftpTree;
  while (<TREE1KG>){
    chomp;
    my @F = split(/\t/,$_);
    
    next if ($_ =~ /evidenceOnly/);
  
    if ($F[0] =~ /^.+\/([^\.]+).+(low_coverage|exome|high_coverage)\.cram$/){
      $URLs{$1}->{$2} = $noPrefix?$F[0]:$ftpPrefix.$F[0];      
    }  
      
    if ($F[0] =~ /phase\S+?\/data\/(\S+?)\/.+\/(.+\.mapped.+(low_coverage|exome|high_coverage)\.\d+\.bam)$/){
      $URLs{$1}->{$3} = $noPrefix?$F[0]:$ftpPrefix.$F[0];      
    }   
  }
}