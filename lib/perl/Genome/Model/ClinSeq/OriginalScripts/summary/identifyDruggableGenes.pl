#!/usr/bin/env genome-perl
#Written by Malachi Griffith

#Given a list of gene symbols in the column X of a TSV file, join drug interaction records for matching gene-to-drug interactions with gene names in column Y and print a merged file

use strict;
use warnings;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use Data::Dumper;

#Required at least one
my $candidates_file = '';
my $interactions_file = '';
my $name_col_1 = '';
my $name_col_2 = '';

GetOptions ('candidates_file=s'=>\$candidates_file, 'interactions_file=s'=>\$interactions_file, 'name_col_1=i'=>\$name_col_1, 'name_col_2=i'=>\$name_col_2);

my $usage=<<INFO;

  Example usage: 
  
  identifyDruggableGenes.pl  --candidates_file=/gscmnt/sata132/techd/mgriffit/hgs/custom/hg1/druggable_genes/drugCandidateGenes.tsv  --name_col_1=1  --interactions_file=/gscmnt/sata132/techd/mgriffit/DruggableGenes/KnownDruggable/DrugBank/query_files/DrugBank_WashU_INTERACTIONS.filtered.3.tsv  --name_col_2=12
  
  Inputs:
  --candidates_file       .tsv file containing candidate gene symbols and any additional files (could be output from mergeDruggableCandidates.pl or something else)
                         Candidates could be mutated genes, CNV gained genes, over-expressed genes
  --name_col_1           Column containing gene symbol for candidates
  --interactions_file    File containing druggable gene interactions.
  --name_col_2           Column containing gene symbol for drug-gene interactions

INFO

unless ($candidates_file && $interactions_file && $name_col_1 && $name_col_2){
  print GREEN, "$usage", RESET;
  exit();
}


#Load the candidate genes - allow for multiple records per candidate
my %candidates;
open (CAN, "$candidates_file") || die "\n\nCould not open candidate file: $candidates_file\n\n";
my $header = 1;
my $header_line1;
my $cc;
while(<CAN>){
  my $record = $_;
  chomp($record);
  if ($header == 1){
    $header = 0;
    $header_line1 = $record;
    next();
  }
  $cc++;
  my @line = split("\t", $record);
  my $gene_id = $line[$name_col_1 - 1];
  if ($candidates{$gene_id}){
    my $record_ref = $candidates{$gene_id}{records};
    $record_ref->{$cc}->{record} = $record;
  }else{
    my %tmp;
    $tmp{$cc}{record} = $record;
    $candidates{$gene_id}{records} = \%tmp;
  }
}
close(CAN);
#print Dumper %candidates;


#Load the candidate genes - allow for multiple records per candidate
my %interactions;
open (INT, "$interactions_file") || die "\n\nCould not open interactions file: $interactions_file\n\n";
$header = 1;
my $header_line2;
my $ic;
while(<INT>){
  my $record = $_;
  chomp($record);
  if ($header == 1){
    $header = 0;
    $header_line2 = $record;
    next();
  }
  $ic++;
  my @line = split("\t", $record);
  my $gene_id = $line[$name_col_2 - 1];
  if ($interactions{$gene_id}){
    my $record_ref = $interactions{$gene_id}{records};
    $record_ref->{$ic}->{record} = $record;
  }else{
    my %tmp;
    $tmp{$ic}{record} = $record;
    $interactions{$gene_id}{records} = \%tmp;
  }
}
close(INT);
#print Dumper %interactions;


#Go through each gene and join matching interactions if present
my $i_count = 0;
print "interaction_count\t$header_line1\t$header_line2\n";
foreach my $gene_id (sort keys %candidates){
  if ($interactions{$gene_id}){
    my $candidates_records = $candidates{$gene_id}{records};
    my $interactions_records = $interactions{$gene_id}{records};

    foreach my $cc (sort {$a <=> $b} keys %{$candidates_records}){
      my $candidate_record = $candidates_records->{$cc}->{record};
      foreach my $ic (sort {$a <=> $b} keys %{$interactions_records}){
        my $interaction_record = $interactions_records->{$ic}->{record};
        $i_count++;
        print "$i_count\t$candidate_record\t$interaction_record\n";
      }
    }
  }
}

exit();
