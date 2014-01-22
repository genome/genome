#!/usr/bin/env genome-perl
# This script will transform a .bed file to exon coordinates
# Written by Malachi Griffith and Rodrigo Goya
use strict;
use warnings;

my $BED_FILE = $ARGV[0] || "-";

my $line_number = 0;
open(BED, "<$BED_FILE") || die("ERROR: could not open junction bed file: '$BED_FILE'\n");
while(<BED>) {	
  chomp; s///;

  #Input format
  #chr3    872615  873354  chr3:872615-873354      23.0    +       872615  873354  0,0,0   2       36,28,  0,711,
  my @f = split(/\t/, $_);
  $line_number++;

  if(!($f[0] =~ m/^chr/)) {next; }
  my $blocks = $f[9];
  my $chr = $f[0];
  my $start = $f[1];
  my $end = $f[2];
  my $name = $f[3];
  my $score = $f[4];
  my $strand = $f[5];
  unless($strand eq "+" || $strand eq "-"){$strand="." };

  my @b_sizes = split(/,/, $f[10]);
  my @b_offsets = split(/,/, $f[11]);
  for(my $b = 0; $b < $blocks; $b++) {
    my $left = $start + $b_offsets[$b-1] + $b_sizes[$b-1];
    my $right = $start + $b_offsets[$b] + 1;
    my $left = $start + $b_offsets[$b] + 1;
    my $right = $start + $b_offsets[$b] + $b_sizes[$b];

    my $code;
    if ($blocks == 1){
      $code = "SingleExon"; #SingleExon
    }elsif($b+1 == 1){
      $code = "MultiExon_First";
    }elsif($b+1 == $blocks){
      $code = "MultiExon_Last";
    }else{
      $code = "MultiExon_Middle";
    }

    #Output format
    #chr12:56412895-56413029 -       meM     CUFF.323719.1   4.7696026179
    #Codes: se (only exon), meUS (first exon of multi-exon transcript), meM (middle exon ...), meDS (last exon ...)
    print "$chr:$left-$right\t$strand\t$code\t$line_number($name)\t$score\n";
  }
}
close(BED);

