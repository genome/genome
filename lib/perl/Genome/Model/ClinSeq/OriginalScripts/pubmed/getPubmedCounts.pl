#!/usr/bin/env genome-perl
#Written by Malachi Griffith

#Take a list of gene names and query pubmed - return the number of abstracts with matches
#Allow an additional string to be combined with the gene name in the search

use strict;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use Data::Dumper;
use XML::Simple;
use above 'Genome'; #remove the 'above' when we turn this into a module

my $gene_list = '';
my $gene_col = '';
my $extra_query = '';
my $sleep = '';
my $db = '';

GetOptions ('gene_list=s'=>\$gene_list, 'gene_col=i'=>\$gene_col, 'extra_query=s'=>\$extra_query, 'sleep=i'=>\$sleep, 'db=s'=>\$db);

unless ($gene_list && $gene_col && $db){
  print RED, "\nRequired input parameter(s) missing or incorrect format!\n\n", RESET;
  print GREEN, "\n\nUsage: getPubmedCounts.pl  --gene_list=mutant_genes.tsv  --gene_col=1  --db=pubmed  --sleep=2", RESET;
  print GREEN, "\n\nUsage: getPubmedCounts.pl  --gene_list=mutant_genes.tsv  --gene_col=1  --db=pubmed  --sleep=2  --extra_query='+AND+Cancer' [Don't forget the plus signs!]", RESET;
  print GREEN, "\n\nGene list is assumed to be tab delimited with gene names in the specified column (1-based)", RESET;
  print GREEN, "\n\nDB must be an actual NCBI DB name (e.g. 'pubmed', 'protein', 'gene', 'omim', etc.)", RESET;
  print GREEN, "\n\tSee full list here: http://eutils.ncbi.nlm.nih.gov/entrez/eutils/einfo.fcgi?", RESET;
  print GREEN, "\n\n'--sleep' is optional and will delay consecutive queries to NCBI\n\n", RESET;
  exit();
}

unless ($sleep){
  $sleep = 0;
}

my %genes;
open(GENES, "$gene_list") || die "\n\nCould not open gene list file: $gene_list\n\n";
my $c = 0;
while(<GENES>){
  $c++;
  chomp($_);
  my @line = split("\t", $_);
  $genes{$line[$gene_col-1]}{line} = $_;
  $genes{$line[$gene_col-1]}{c} = $c;
}
close(GENES);

foreach my $gene (sort {$genes{$a}{c} <=> $genes{$b}{c}} keys %genes){

  my $term = $gene;

  if ($extra_query){
    $term .= "$extra_query";
  }

  my $xml_file = "temp.xml";
  #print BLUE, "\n\nGene: $gene", RESET;
  my $wget_cmd = "wget \"http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db="."$db"."&term=$term\" -O $xml_file 2>/dev/null";
  #print YELLOW, "\n$wget_cmd", RESET;
  Genome::Sys->shellcmd(cmd => $wget_cmd);

  my $xs1 = XML::Simple->new();
  my $doc = $xs1->XMLin($xml_file);

  my $abstract_count = $doc->{'Count'};
  #print BLUE, "\n\tFound $abstract_count abstracts with the term: $term", RESET;

  print "$genes{$gene}{line}\t$abstract_count\n";
  sleep $sleep;
}

#print "\n\n";

exit();



