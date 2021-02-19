use strict;

use Getopt::Long;
use List::MoreUtils qw(uniq);

GetOptions ('zf=s'  => \(my $fileZF),
	          'all=s' => \(my $fileAllele));

################################################################################
my %zfs;
open IN, $fileAllele;
while(<IN>){
  chomp;
  my ($name,$zf) = split(/\t/,$_);
  $zfs{$name} = $zf;
}
close IN;

################################################################################
my %zfseq;
open IN, $fileZF;
while(<IN>){
  chomp;
  my ($name,$zf) = split(/\t/,$_);
  $zfseq{$name} = $zf;
}
close IN;

################################################################################
my (%editdist_scores,%alignment_scores);

for my $zfa_i(sort keys (%zfs)){

  ($editdist_scores{$zfa_i.":".$zfa_i},$alignment_scores{$zfa_i.":".$zfa_i}) = getPRDM9alignmentScore($zfa_i,$zfa_i) unless($alignment_scores{$zfa_i.":".$zfa_i});

  for my $zfa_j(sort keys (%zfs)){
    my $combo    = $zfa_i.":".$zfa_j;
    my $revcombo =  $zfa_j.":".$zfa_i;

    if ($alignment_scores{$revcombo}){
      ($editdist_scores{$combo},$alignment_scores{$combo}) = ($editdist_scores{$revcombo},$alignment_scores{$revcombo});
    }

    if (!$alignment_scores{$combo}){
      ($editdist_scores{$combo},$alignment_scores{$combo}) = getPRDM9alignmentScore($zfa_i,$zfa_j);
    }

    print join("\t",$zfa_i,
                    $zfa_j,
                    $editdist_scores{$combo},
                    $alignment_scores{$combo},
                    sprintf("%4.3f",$alignment_scores{$combo}/$alignment_scores{$zfa_i.":".$zfa_i}*100))."\n";
  }
}

################################################################################
sub getPRDM9alignmentScore{
  my ($zfa_allele_1, $zfa_allele_2) = @_;

  ## Get all ZFs for the two alleles (alternative to split coz of 2 chars)
  my @zf1 = unpack '(a2)*', $zfs{$zfa_allele_1};
  my @zf2 = unpack '(a2)*', $zfs{$zfa_allele_2};

  ## Get single letter codes for each ZF
  my $zf_cnt = 0;
  my %zf_codes;
  my %sub_matrix;

  for my $uZF(uniq(sort (@zf1,@zf2))){
    my $newZFcode = chr(65+($zf_cnt++));
    $zf_codes{$uZF} = $newZFcode;
    $zfseq{$newZFcode} = $zfseq{$uZF};
  }

  ## Convert to new single letter code
  for my $i(0..$#zf1){$zf1[$i] = $zf_codes{$zf1[$i]}}
  for my $i(0..$#zf2){$zf2[$i] = $zf_codes{$zf2[$i]}}

  my $tmp_folder = './';
  my $tmp_stem   = $tmp_folder.int(rand()*100000000000000);
 
  my $matrix_file = $tmp_stem.'_SMAT';

  my $zf1FA       = $tmp_stem.'_zfa1.fa';
  my $zf2FA       = $tmp_stem.'_zfa2.fa';

  make_substitution_matrix($zf_cnt-1,$matrix_file,\%sub_matrix);
  make_FASTA($zfa_allele_1,join("",@zf1),$zf1FA);
  make_FASTA($zfa_allele_2,join("",@zf2),$zf2FA);

  my @aln   = `needle -asequence $zf1FA -bsequence $zf2FA -gapopen 1 -gapextend 0 -datafile $matrix_file -aformat3 fasta -stdout -auto`;
  my $score = `needle -asequence $zf1FA -bsequence $zf2FA -gapopen 1 -gapextend 0 -datafile $matrix_file -aformat3 score -stdout -auto |head -n1`;
  $score  =~ s/^.+\(([\d\.]+)\).*$/$1/;

  for my $i(0..3){chomp $aln[$i]};

  my ($edit_distance,$aln_score) = ndiffSM($aln[1],$aln[3],$score,\%sub_matrix);

  system("rm $tmp_stem*");
  return($edit_distance,$aln_score);
}

################################################################################
sub make_substitution_matrix{
  my ($nZF,$mout,$sub_matrix) = @_;

  open SM, '>', $mout;

  my $header = " ";
  for my $i(0..$nZF){$header .= padMe(chr(65+$i))}
  print SM $header."\n";

  for my $i(0..$nZF){
    my $smLine = chr(65+$i);

    for my $j(0..$nZF){
      my $dist = 1;
      if ($i ne $j){
        $dist = ndiff($zfseq{chr(65+$i)},$zfseq{chr(65+$j)});
        $dist = ($dist > 0)?$dist*-1:1;
      }

      $smLine .= padMe($dist);

      ## Keep in RAM too
      $$sub_matrix{chr(65+$i).chr(65+$j)} = $dist;
    }
    print SM $smLine."\n";
  }
  close SM;
}

################################################################################
sub make_FASTA{
  my ($name, $seq, $fa) = @_;

  open ZF, '>', $fa;
  print ZF ">$name\n$seq";
  close ZF;
}

################################################################################
sub ndiff{
  my ($str1,$str2) = @_;

  my $mask = $str1 ^ $str2;
  my $diff_cnt;
  $diff_cnt++ while ($mask =~ /[^\0]/g);

  return($diff_cnt);
}

################################################################################
sub ndiffSM{
  my ($str1,$str2,$aln_score,$sub_matrix) = @_;

  my $mask = $str1 ^ $str2;

  my $score;

  $score -= 1 while ($str1 =~ m/\-+/g);
  $score -= 1 while ($str2 =~ m/\-+/g);

  while ($mask =~ /[^\0]/g) {
      my $zf1 = substr($str1,$-[0],1);
      my $zf2 = substr($str2,$-[0],1);

      next if ($zf1 eq "-" || $zf2 eq "-");

      $score += $$sub_matrix{$zf1.$zf2};

      #print join("\t",$zf1,$zf2,$$sub_matrix{$zf1.$zf2})."\n";
  }

#  print join("\t",$str1,$str2,$score*-1,$aln_score      )."\n";
#  print "--------------------------------\n";
#  print join("\t",$str1,$score*-1)."\n";
#  print join("\t",$str2,$score*-1)."\n";

  chomp $score;
  chomp $aln_score;

  return($score?$score:0, $aln_score);
}

################################################################################
sub padMe{
  my ($val) = shift;

  while (length($val) < 4){$val = " ".$val}
  return($val);
}
