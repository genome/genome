#!/usr/bin/env genome-perl -w
#Written by Malachi Griffith

use strict;
use Getopt::Long;
use Data::Dumper;
use Term::ANSIColor qw(:constants);

#Start with nanostring data - one line per refseq
my $nanostring_data = "/gscmnt/sata132/techd/twylie/Stratagene_Colon_RNA-seq/cuffDiffCompare/Nanostring/20110504-Colon-Cancer-Reference.txt";
my %nano;
open (NANO, $nanostring_data) || die "\n\nCould not open nanostring data\n\n";
while(<NANO>){
  chomp($_);
  my @line = split("\t", $_);
  if ($_ =~ /Endogenous/){
    my $tmp = $line[2];
    if ($tmp =~ /(\w+)\.\d+/){
      my $id = $1;
      if ($line[14] =~ /(\d+)/){
        $line[14] = $1;
        $nano{$id}{line} = \@line;
      }
    }
  }
}
close(NANO);
#print Dumper %nano;


#Get the Cufflinks data - one line per UCSC id
my $cufflinks_data = "/gscmnt/sata132/techd/twylie/Stratagene_Colon_RNA-seq/cuffDiffCompare/4_x_4_trials/genes.fpkm_tracking.simple.txt";
my %cuff;
open (CUFF, $cufflinks_data) || die "\n\nCould not open cufflinks data\n\n";
while(<CUFF>){
  chomp($_);
  my @line = split("\t", $_);
  if ($line[2] =~ /\d+/){
    $cuff{$line[0]}{line} = \@line;
  }
}
close(CUFF);
#print Dumper %cuff;

#Get the refseq to ucsc mappings
my $mapping_data = "/gscmnt/sata132/techd/twylie/Stratagene_Colon_RNA-seq/cuffDiffCompare/Nanostring/UCSC_vs_Refseq_IDs.tsv";
my %map;
open (MAP, $mapping_data) || die "\n\nCould not open map data\n\n";
while(<MAP>){
  chomp($_);
  my @line = split("\t", $_);
  my $refseq = $line[2];
  my $ucsc = $line[0];
  if ($refseq){
    if ($map{$refseq}){
      my $ucsc_ref = $map{$refseq};
      $ucsc_ref->{$ucsc}=1;
    }else{
      my %tmp;
      $tmp{$ucsc}=1;
      $map{$refseq} = \%tmp;
    }
  }
}
close(MAP);
#print Dumper %map;

#Print out a new nanostring data file mapped to the cuff
#Since one refseq is often associated with multiple UCSC IDs, duplicate the Nanostring data where neccessary
foreach my $refseq (sort keys %nano){
  if ($map{$refseq}){
    my $ucsc_ids_ref = $map{$refseq};
    foreach my $ucsc (sort keys %{$ucsc_ids_ref}){
      my @nl = @{$nano{$refseq}{line}};
      if ($cuff{$ucsc}){
        my @cuff = @{$cuff{$ucsc}{line}};
        print "$nl[2]\t$ucsc\t$nl[3]\t$nl[4]\t$nl[5]\t$nl[6]\t$nl[7]\t$nl[8]\t$nl[9]\t$nl[10]\t$nl[11]\t$nl[12]\t$nl[13]\t$nl[14]\t$cuff[1]\t$cuff[2]\t$cuff[3]\t$cuff[4]\n";
      }
    }


  }
}









