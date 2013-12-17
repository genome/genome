#!/usr/bin/env genome-perl
# Plot breakdancer SVs using pairoscope and Xian's copy number graph tool
use strict;
use warnings;
use Getopt::Std;

my %opts = (
	    l=>5000,
	    q=>1
	   );
getopts('l:ACa:b:c:d:D:2q:',\%opts);
die("
Add annotation (UCSC gene, dbSNP and others) to a text file
Usage:   pairoscope.pl <a BreakDancer output file>\n
Options:
         -2       Plot gmt copy-number graph (instead of pairoscope), run on 64bit machine
         -l INT   Proximity of breakpoints to gene annotations [$opts{l}] bp
         -D STR   Output directory
         -A       Annotate Merged SV assembly file
         -C       Plot Captured bams
         -q       Minimal Mapping Quality [$opts{q}]
         -a INT   Which column contains the chromosome of the breakpoint 1 (1-based)
         -b INT   Which column contains the position of the  breakpoint1 (1-based)
         -c INT   Which column contains the chromosome of the breakpoint 2 (1-based)
         -d INT   Which column contains the position of the breakpoint 2 (1-based)
\n") unless (@ARGV);

if($opts{A}){
  $opts{a}=2; $opts{b}=4; $opts{c}=5; $opts{d}=6;  #SV assembly format
}
$opts{D}='./' if(!defined $opts{D});

my $header;
foreach my $fin(@ARGV){
  my $FIN;
  if($fin eq '-'){
    $FIN=\*STDIN;
  }
  else{
    open(IN,"<$fin") || die "unable to open $fin\n";
    $FIN=\*IN;
  }
  while(<$FIN>){
    chomp;
    if(/^\#/){
      $header=$_ if(/^\#/ && /Chr/i);
      next;
    }
    chomp;
    my @u=split /\s+/;
    my ($chr,$start,$ori1,$chr2,$end,$ori2,$type,$size,$score,$nreads,$nread_lib,$AF,@extra)=@u;
    next if(defined $opts{o} && $chr ne $opts{o});

    $chr=$u[$opts{a}-1] if(defined $opts{a});
    $start=$u[$opts{b}-1] if(defined $opts{b});
    $start=~s/\D.*//gi;
    $chr2=$u[$opts{c}-1] if(defined $opts{c});
    $end=$u[$opts{d}-1] if(defined $opts{d});
    $end=~s/\D.*//gi;
    next unless($start=~/^\d+/ && $end=~/^\d+/ && defined $chr);
    if($opts{2}){  #copy-number graph
      my $fname=join('.',$chr,$start,$chr2,$end,'copy-number');
      my $tbam=($opts{C})?"tumor.capture.bam":"tumor.bam";
      my $nbam=($opts{C})?"normal.capture.bam":"normal.bam";
      my $cmd=sprintf "gmt copy-number graph --chromosome=%s --start=%d --end=%d --output-dir=./ --name=%s --tumor-bam-file=$tbam --normal-bam-file=$nbam",$chr,$start,$end,$fname;
      print "$cmd\n";
      `bsub -N -u plotSV\@$ENV{GENOME_EMAIL_DOMAIN} -M 8000000 -R "select[type=LINUX64] select[mem>8000] rusage[mem=8000]" -J COPYNUMBER -oo $fname.log $cmd`;
    }
    else{  #pairoscope
      my @bams=($opts{C})?('tumor.capture','normal.capture'):('tumor','normal');
      foreach my $s(@bams){
	my $fname=$opts{D} . join('.',$chr,$start,$chr2,$end,$s,'pairoscope.png');
	my $cmd=sprintf "~dlarson/src/pairoscope/trunk/pairoscope -q $opts{q} -g /gscmnt/sata831/info/medseq/dlarson/annotation_sam_files/new_annotation_sorted.bam -o %s %s.bam %s %d %d %s.bam %s %d %d", $fname,$s,$chr,$start-$opts{l},$start+$opts{l},$s,$chr2,$end-$opts{l},$end+$opts{l};
        print "$cmd\n";
	system($cmd);
      }
    }
  }
}
