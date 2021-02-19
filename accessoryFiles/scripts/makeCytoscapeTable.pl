use strict;

my %alleleOK;
my %potentialParents;

my $includeAll = $ARGV[0];

if (-e 'prdm9_haplotypes.final.tab'){
    open IN, 'prdm9_haplotypes.final.tab';

    while (<IN>){
        chomp;
        $_ =~ s/Av:s:0053:M1S:A-A/N/g;
        $_ =~ s/(\S+v):[cs]:(\d+)\S+/$1:$2/g;
        $_ =~ s/baudat_(\S+)/b$1/g;

        my @F = split(/\t/,$_);

        $alleleOK{$F[4]}++ if ($F[4] !~ /(Unk|allele|Paren|PrZFA|noZFA|NA)/);
        $alleleOK{$F[5]}++ if ($F[5] !~ /(Unk|allele|Paren|PrZFA|noZFA|NA)/);
    }
    close IN;
}

open REC, 'PrZFA_recombinants_ALL.tab';

print join("\t","PrZFA","parent","7-N")."\n";
while (<REC>){
    chomp;
    $_ =~ s/Av:s:0053:M1S:A-A/N/g;
    $_ =~ s/(\S+v):[cs]:(\d+)\S+/$1:$2/g;
    $_ =~ s/baudat_(\S+)/b$1/g;

    my @F = split(/\t/,$_);

    if ($includeAll || $alleleOK{$F[0]} || $F[0] =~ /^([A-Z]|[LM]\d+|b[A-Z]|pt[A-Z])$/i){
        if ($includeAll || $alleleOK{$F[1]} || $F[1] =~ /^([A-Z]|[LM]\d+|b[A-Z]|pt[A-Z])$/){

            $alleleOK{$F[0]}++;
            $alleleOK{$F[1]}++;
            if ($_ =~ /END/){
                print join("\t",@F[0..1],7-$F[3])."\n";
                $potentialParents{$F[0]}++;
            }
        }
    }
}
close REC;

for my $allele(sort keys(%alleleOK)){
    print join("\t",$allele,$allele,0)."\n" unless ($potentialParents{$allele});
}
