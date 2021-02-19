use strict;

my ($idFile, $type, $fqPath) = @ARGV;

my (  %reads, %ids);
open IN, $idFile;
while (<IN>){
  chomp;
  next if ($_ =~ /bc11.+individu/ || $_ =~ /unk/i);
  my ($id1,$id2,$bc11,$bc12,$bc21,$bc22,$name,$individual,$pop,$popid,$sex) = split(/\t/,$_);
  my $name = $individual."_".$pop;
  $reads{$id1} = $name;
  $reads{$id2} = $name;
  $ids{$name}++;
}
close IN;

my %fh;
for my $name(sort(keys %ids)){
  open ($fh{$name}, '>>', $name.".$ARGV[1].fastq");
}

opendir(DIR, $fqPath);

while (my $fl = readdir(DIR)){
    next unless ($fl =~ /(fastq|fq)/);
    my $FQ;
    if ($fl =~ /\.gz$/){
      open $FQ, '-|', 'zcat '.$ARGV[2].'/'.$fl;
    }else{
      open $FQ, '-|', 'cat '.$ARGV[2].'/'.$fl;
    }

    my @a;

    ## Read FASTQ file 4 lines at a time
    while ( (@a) = map { (scalar(<$FQ>) or ()) } (0..3) ) {
      $a[0] =~ /^\@(\S+)/;
      my $read = $1;

      if ($reads{$read}){
        my $name = $reads{$read};
        $a[0] =~ s/^(\@\S+)\s+/$1 individual=$name basecalls=guppy3.1.5 basecalltype=$type /;
        my $FH = $fh{$name};
        print $FH (@a);
      }
    }
  }
