#!/usr/bin/env genome-perl
# Malachi Griffith and Rodrigo Goya

# This script will take a BED stream and output the a list of junctions
# Each unique junction is defined by chr:start-end(strand)
# The number of times each is observed will be noted
# Example commands:
# cat CCDS.Genes.bed | combineKnownJunctionBeds.pl > CCDS.Genes.junc
# cat CCDS.Genes.bed Ensembl.Genes.bed MGC.Genes.bed Refseq.Genes.bed UCSC.Genes.bed Vega.Genes.bed | combineKnownJunctionBeds.pl > ALL_COMBINED.junc

#Input format:
#Standard BEDs downloaded from UCSC genome browser

#Output format:
#chr:start-end(strand)  chr start end strand  count transcript_id_list gene_id_list

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use File::Basename;
use IO::File;
use Benchmark;

my $map_file = '';

GetOptions ('map_file=s'=>\$map_file);
if ($map_file){
  #parameter ok
}else{
  print RED, "\nRequired input parameter(s) missing\n\n", RESET;
  print GREEN, "e.g. cat CCDS.Genes.bed | combineKnownJunctionBeds.pl  --map_file=transcript_to_gene/CCDS.Genes.map > CCDS.Genes.junc\n\n", RESET;
  print GREEN, "e.g. cat CCDS.Genes.bed Ensembl.Genes.bed MGC.Genes.bed Refseq.Genes.bed UCSC.Genes.bed Vega.Genes.bed | combineKnownJunctionBeds.pl --map_file=transcript_to_gene/ALL.Genes.map > ALL_COMBINED.junc\n\n", RESET;
  exit();
}

#Use transcript to gene mapping files
#This file contains transcript ID to gene mapping info
#my $map_file = "/projects/alexa2/hmmSplicer/ReferenceJunctions/hg18/transcript_to_gene/ALL.Genes.map";
my %transcripts;
open(MAP, "$map_file") || die "\n\nCould not open map file: $map_file\n\n";
while(<MAP>){
  chomp($_);
  my @line = split("\t", $_);
  my $trans_id = $line[0];
  my $gene_id_string = $line[1];

  if ($gene_id_string eq "n/a"){
    next();
  }
  my @gene_ids = split(",", $gene_id_string);
  foreach my $gid (@gene_ids){
    if ($gid eq "n/a"){
      next();
    }

    if ($transcripts{$trans_id}){
      my $gids_ref = $transcripts{$trans_id}{gids};
      $gids_ref->{$gid} = 1;
    }else{
      my %tmp;
      $tmp{$gid} = 1;
      $transcripts{$trans_id}{gids} = \%tmp;
    }
  }
}
close(MAP);

foreach my $tid (keys %transcripts){
  my $gids_ref = $transcripts{$tid}{gids};
  my $gid_string = '';
  my $first = 1;
  foreach my $gid (keys %{$gids_ref}){
    if ($first){
      $gid_string .= "$gid";
      $first = 0;
    }else{
      $gid_string .= ",$gid";
    }
  }
  $transcripts{$tid}{gid_string} = $gid_string;
}

#Parse junctions out of a BED stream
my $BED_FILE = $ARGV[0] || "-";
my $line_number = 0;
open(BED, "<$BED_FILE") || die("ERROR: could not open junction bed file: '$BED_FILE'\n");
my %junc;
while(<BED>) {
  chomp; s///;
  #chr3	189881601	189888428	JUNC00064306	4	+	189881601	189888428	255,0,0	2	20,29	0,6798
  my @f = split(/\t/, $_);
  $line_number++;

  if(!($f[0] =~ m/^chr/)) {next; }
  my $blocks = $f[9];
  if($blocks == 1) { next; }

  my $chr = $f[0];
  my $start = $f[1];
  my $trans_id = $f[3];
  my @b_sizes = split(/,/, $f[10]);
  my @b_offsets = split(/,/, $f[11]);
  my $strand = $f[5];

  for(my $b = 1; $b < $blocks; $b++) {
    my $left = $start + $b_offsets[$b-1] + $b_sizes[$b-1];
    my $right = $start + $b_offsets[$b] + 1;
    my $junc_name = "$chr:$left-$right";
    my $j_id = "$chr:$left-$right($strand)";
    
    if ($junc{$j_id}){
      $junc{$j_id}{count}++;
      $junc{$j_id}{trans_ids} .= ",$trans_id";
    }else{
      $junc{$j_id}{name} = $junc_name;
      $junc{$j_id}{count} = 1;
      $junc{$j_id}{chr} = $chr;
      $junc{$j_id}{left} = $left;
      $junc{$j_id}{right} = $right;
      $junc{$j_id}{strand} = $strand;
      $junc{$j_id}{trans_ids} = "$trans_id";
    }
  }
}
close(BED);

#Determine the gene ids associated with the list of transcript ids
foreach my $jid (sort {$a cmp $b} keys %junc){
  my $trans_id_string = $junc{$jid}{trans_ids};
  my @trans_ids = split(",", $trans_id_string);
  my %gene_ids;
  foreach my $tid (@trans_ids){
    my $gid_string = $transcripts{$tid}{gid_string};
    my @gene_ids = split(",", $gid_string);
    foreach my $gid (@gene_ids){
      $gene_ids{$gid}=1;
    }
  }
  my $final_gid_string = '';
  my $first = 1;
  foreach my $gid (sort keys %gene_ids){
    if ($first){
      $final_gid_string .= "$gid";
      $first = 0;
    }else{
      $final_gid_string .= ",$gid";
    }
  }
  unless($final_gid_string =~ /\S+/){
    $final_gid_string = "n/a";
  }
  $junc{$jid}{final_gid_string} = $final_gid_string;
}

foreach my $jid (sort {$a cmp $b} keys %junc){
  print "$jid\t$junc{$jid}{chr}\t$junc{$jid}{left}\t$junc{$jid}{right}\t$junc{$jid}{strand}\t$junc{$jid}{count}\t$junc{$jid}{final_gid_string}\t$junc{$jid}{trans_ids}\n";
}


