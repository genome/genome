#!/usr/bin/env genome-perl
# This script generates read counts and copy number from a single bam file,
# which can be used for copy number segmentation

use strict;
use warnings;
use Getopt::Std;
use Statistics::Descriptive;

my $version="BAM2CN-0.0.1r2";
my %opts = (w=>10000, q=>0, n=>0.25, r=>1);
my %opts1;
getopts('w:q:c:n:r:p', \%opts1);
die("
Usage:   BAM2CN.pl <a bam file>
Output read depth based copy number in a single genome in one or more map/bam files
Options:
         -c STRING   Chromosome number
         -n FLOAT    ratio diverged from median, used to find copy number neutral region
         -w INT      Window size [$opts{w} bp]
         -q INT      MAQ mapping quality cutoff [$opts{q}]
         -r FLOAT    downsampling ratio [$opts{r}]
         -p          Only examine properly mapped reads
Version: $version\n
") unless (@ARGV);

my $options='';
foreach my $opt(keys %opts1){
  $options.=$opt.$opts1{$opt};
  $opts{$opt}=$opts1{$opt};
}

######################## Read In Configuration file #######################
my $fbam=$ARGV[0];

######################## Compute read counts in sliding windows ########################
my %data;
my %Coor;
my %Chrs;

my $cmd="samtools view $fbam";
$cmd.=" $opts{c}" if(defined $opts{c});
open(MAP,"$cmd |") || die "unable to open $fbam\n";
my ($ppos,$nread_buf);
my $pchr=0;
my $idx=-1;
while(<MAP>){
  chomp;
  my $rr=rand();
  next if($rr>$opts{r});  #down sampling reads
  my @t=split;
  my ($chr,$pos,$mapqual,$dist,$flag);
  $flag=0;
  $mapqual=$t[4];
  $chr=$t[2];
  $pos=$t[3];
  $flag=1 if(!defined $opts{p} || ($t[1]=~/P/ || $t[1] & 0x0002));

  next unless ($mapqual=~/^\d+$/ && $mapqual>=$opts{'q'} && $flag && $chr=~/\S+/ && $pos=~/^\d+$/);
  next if(defined $opts{'c'} && $opts{'c'} ne $chr);

  if($chr ne $pchr){  # a new chromosome
    # reset
    $Chrs{$chr}=1;
    $pchr=$chr;
    $nread_buf=0;
    $ppos=0;
    $idx=-1;
  }
  $nread_buf++;
  if($pos>$ppos+$opts{'w'}){
    do{
      $ppos+=$opts{'w'};
    } until($pos<$ppos+$opts{'w'});

    #register
    $idx++;
    ${$data{$chr}}[$idx]+=$nread_buf;
    ${$Coor{$chr}}[$idx]=$ppos;
    $nread_buf=0;
  }
}
close(MAP);


my @chrs=sort keys %Chrs;
die "No reads passed QC.\n" if($#chrs<0);

my $median=Statistics::Descriptive::Full->new();
foreach my $chr(@chrs){
  my $tmpchr=$chr; $tmpchr=~s/chr//i;
  next unless($tmpchr=~/^\d+/);  #skip non-autosome
  my $chr_median=&Get_Median($data{$chr});
  #print "#Chr${chr}_Median:$chr_median\n";
  $median->add_data($chr_median);
}
my $md=$median->median();

print "#WholeGenome_Median:$md\n";
my $num_CN_neutral_pos=0;
my $NReads_CN_neutral=0;
foreach my $chr(@chrs){
  next unless (defined $data{$chr});
  my $Nwi=$#{$data{$chr}};

  for(my $i=0;$i<=$Nwi;$i++){
    my $f2x=1;
    next unless (defined ${$data{$chr}}[$i]);

    $f2x=0 if(${$data{$chr}}[$i]<$md*(1-$opts{n}) || ${$data{$chr}}[$i]>$md*(1+$opts{n}));
    next if(! $f2x);
    $num_CN_neutral_pos++;
    $NReads_CN_neutral+=${$data{$chr}}[$i];
  }
}

#subtract the normal from the tumor
my $depth2x=($num_CN_neutral_pos>10)?$NReads_CN_neutral/$num_CN_neutral_pos:$md;
printf "#2xReadCount:%.2f in %d bp\n",$depth2x,$opts{w};
print "CHR\tPOS\tReadCount\tCopyNumber\n";

foreach my $chr(@chrs){
  next unless (defined $data{$chr});
  my $Nwi=$#{$data{$chr}};

  for(my $i=0;$i<=$Nwi;$i++){
    next unless (defined ${$data{$chr}}[$i]);
    my $cn=${$data{$chr}}[$i]*2/$depth2x;
    my $poschr=${$Coor{$chr}}[$i];
    printf "%s\t%d\t%.6f\t%.6f\n",$chr,$poschr,${$data{$chr}}[$i],${$data{$chr}}[$i]*2/$depth2x;
  }
}

sub Get_Median {
  my $rpole = shift;
  my @pole = @$rpole;
  my $ret;

  @pole=sort {$a<=>$b} @pole;
  if( (@pole % 2) == 1 ) {
    $ret = $pole[((@pole+1) / 2)-1];
  } else {
    $ret = ($pole[(@pole / 2)-1] + $pole[@pole / 2]) / 2;
  }
  return $ret;
}
