use strict;
use Bio::SeqIO;

system(`esearch -db nuccore -query "Prdm9 AND Mus musculus[Organism]" |efetch -format fasta >mmPRDM9.fa`);

my $fa  = Bio::SeqIO->new(-file => "mmPRDM9.fa" , -format => 'fasta');
#my $fa  = Bio::SeqIO->new(-file => "allfa.txt" ,
#                          -format => 'fasta');

my %zfs ;
my $nzf = 64;
my $alleleCnt = 0;

open TAB, '>', 'mousePrZFA.tab';
open FA,  '>', 'mousePrZFA.fa';
while ( my $seq = $fa->next_seq() ) {
   if ($seq->desc =~ /prdm9|histone\-lysine\sN\-methyltransferase|hybrid\ssterility\s1\sregion/i){
      my $name = $seq->desc;

      $name =~ s/(,|PREDICTED|strain|and|partial|cds|isolate|PR\-domain|protein|PR\sdomain|containing|complete|gene|exon|musculus|Mus|:)//gi;
      $name =~ s/^\s+//;
      $name =~ s/\s+(\S)/_$1/gi;
      $name =~ s/for_12_sequence_//gi;

      my $s = $seq->seq;
      if ($s =~ s/^.+?(TGTG.{80}TGTG.{80}TGTG.+)$/$1/){
         if ($s =~ s/^(.+TGTG.{80}TGTG.{80}).+$/$1/){
            if ($name !~ /target/){
               my $out;
               my $alleleZFs;
               my $allele = ++$alleleCnt;
               my $strain = getStrain($seq->desc);
               while ($s =~ s/^(.{84})//){
                  my $zf = $1;
                  $out .= $zf."\n";
                  print TAB $strain."_".$seq->display_id."\t$zf\n";
               }
               print FA "\>".$strain."_".$seq->display_id."\n".$out ;
            }
         }
      }
   };
}

sub getStrain {
   my $txt = shift;
   return "MOL" if ($txt =~ /(molossinus|MOLF\/EiJ)/i);
   return "DOM" if ($txt =~ /(domesticus|C57b)/i);
   return "SPR" if ($txt =~ /spretus/i);
   return "CST" if ($txt =~ /castaneus/i);
   return "MUS" if ($txt =~ /musculus\smusculus/i);
   return "UNK" ;
}
