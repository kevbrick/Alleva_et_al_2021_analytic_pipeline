#!/usr/bin/perl

use strict; 

my ($bam, $chunk_size) = @ARGV;

my $stem = $bam; $stem =~ s/.bam/_subset/;

system("samtools view -H $bam >head.sam");

my $cmd = "samtools view $bam";

open my $PIPE, '-|', $cmd; 

my $bamcount = 0;

my $sam_out = $stem."_".($bamcount++).".sam";
system("samtools view -H $bam >$sam_out");
open SAM, '>>', $sam_out;

my $n_current = 0;
my %ZMWs;

while(<$PIPE>){
	chomp; 
	$_ =~ /^.+\szm\S+:(\S+)\s.+$/;
	my $ZMW = $1; 
	$n_current++ unless ($ZMWs{$ZMW}++);
	print SAM $_."\n";
	if ($n_current > $chunk_size){
		close SAM; 
		my $bam = $sam_out; $bam =~ s/.sam/.bam/;
		system("samtools view -Shb $sam_out >$bam");
		system("rm $sam_out");

		$sam_out = $stem."_".($bamcount++).".sam";
		
		system("samtools view -H $bam >$sam_out");
		open SAM, '>>', $sam_out;
		$n_current = 0;
	}
}

close SAM; 
my $bam = $sam_out; $bam =~ s/.sam/.bam/;
system("samtools view -Shb $sam_out >$bam");
system("rm $sam_out");
