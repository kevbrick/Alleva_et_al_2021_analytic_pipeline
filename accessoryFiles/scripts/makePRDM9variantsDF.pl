use strict; 

## -----------------------------------------------------------------------------
open IN, 'humanPRDM9alleles.BloodandSpermVariants.txt';

my (%keep, %done, %mangt); 

## Get all variants per individual from Jeffreys 2013
while (<IN>){
  chomp; 
  my @F = split(/\t/,$_);
  #M9S:C-L12	:A:B:C:D:D:C:C:F:K:H:L:H:H:I:J	Cv:s:0208:M7S:C-C
  $F[0] =~ /^(M\d+)([S|B]):(\S+)\-(\S+)/;
  my ($manOnly, $man, $p1, $p2) = ($1, $1.$2, $3, $4);
  
  $mangt{$manOnly} = "$p1\/$p2";
  
  $keep{join("\|",$F[2],$p1,$p1)} .= $man."|" unless ($done{join("\|",$F[2],$p1,$p1,$man)}++);
  $keep{join("\|",$F[2],$p2,$p2)} .= $man."|" unless ($done{join("\|",$F[2],$p2,$p2,$man)}++);
  $keep{join("\|",$F[2],$p1,$p2)} .= $man."|" unless ($done{join("\|",$F[2],$p1,$p2,$man)}++);
  $keep{join("\|",$F[2],$p2,$p1)} .= $man."|" unless ($done{join("\|",$F[2],$p2,$p1,$man)}++);
}
close IN; 

## -----------------------------------------------------------------------------
open AC, 'PrZFA_alleles.details.txt';

my (%ACtype); 
my (%inpop); 

## Get all A-type/C-type designations
while (<AC>){
  chomp; 
  my ($id,$shortid,$zf,$pop,$new,$nt,$aa,$aaBS,$nToA,$nToC,$AC) = split(/\t/,$_);
  
  
  $ACtype{$id}      = $AC;
  $ACtype{$shortid} = $AC;
  
  $inpop{$id}++ if ($id =~ /^ABCDEFGHIL(\d*)$/);
}
close AC; 

## -----------------------------------------------------------------------------
open POP, 'prdm9_haplotypes.final.tab';

## Get all in-population alleles
while (<POP>){
  chomp; 
  my @F = split(/\t/,$_);
  
  $inpop{$F[4]}++;
  if ($F[4] =~ /^(.+?):\S+(:\d+):/){;
    my $shortID = $1.$2;
    $inpop{$shortID}++;
  }
  
  $inpop{$F[5]}++;
  if ($F[5] =~ /^(.+?):\S+(:\d+):/){
    my $shortID = $1.$2;
    $inpop{$shortID}++;
  }
}
close POP; 

## -----------------------------------------------------------------------------
open my $IN,  '-|', 'cat PrZFA_relatedness*tab';
open TJEF,    '>', 'parentage_Jeffreys_2013_all.txt';
open TPOP,    '>', 'parentage_Pop_alleles.txt';
open TPOPALL, '>', 'parentage_Pop_alleles_ONLY.txt';

print TJEF    join("\t","id","p1","p2","type","code","switchkey","nswitches","manid","sb","mangt")."\n";
print TPOP    join("\t","id","p1","p2","type","code","switchkey","nswitches","childAC","p1AC","p2AC")."\n";
print TPOPALL join("\t","id","p1","p2","type","code","switchkey","nswitches","childAC","p1AC","p2AC")."\n";

while (<$IN>){
  chomp; 
  my @F = split(/\t/,$_);
  my $type = join("\|",@F[0..2]);
  next if ($_ =~ /^ID/);
  if ($inpop{$F[0]}){
    print TPOP join("\t",@F[0..6],$ACtype{$F[0]},$ACtype{$F[1]},$ACtype{$F[2]})."\n" ;
    if ($inpop{$F[1]} && $inpop{$F[2]}){
      print TPOPALL join("\t",@F[0..6],$ACtype{$F[0]},$ACtype{$F[1]},$ACtype{$F[2]})."\n" ;
    }
  }
  
  for my $man ($keep{$type} =~ /(M\d+[SB])\|/g){
    $man =~ /M(\d+)([SB])/;
    my ($manid, $num, $sb) = ("M$1",$1,$2);
    $sb = ($sb eq "S"?"Sperm":"Blood");
    my @out = @F; 
    push @out, $manid; 
    push @out, $sb; 
    push @out, $mangt{$manid}; 
    $out[6] = 99 if ($_ =~ /NORECOMBINANT/);
    print TJEF join("\t",@out)."\n";
  }

}
close IN;
close TJEF;
close TPOP;
close TPOPALL;