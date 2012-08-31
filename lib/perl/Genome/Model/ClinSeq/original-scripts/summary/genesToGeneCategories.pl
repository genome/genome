#!/usr/bin/env genome-perl
#Written by Malachi Griffith

#This script takes an input gene list and annotates it against various lists of gene setting values as 0/1
#e.g. I have a list of 1000 genes, and I want to know which is a kinase, transcription factor, etc.

#Input parameters / options
#Input file (containing gene names)
#Gene name column.  Column number containing gene symbols
#Symbol lists to annotate with (display a list of gene symbol lists to select from and the location being queried)
#Output file

#1.) Take an input file with gene names in it
#2.) Get the gene name column from input
#3.) 'fix' gene names to Entrez official gene symbols
#4.) Load the symbols lists (fixing gene names on each)
#5.) Intersect the gene names in the input list with each symbol list
#6.) Print output file with annotations and new column headers appended

use strict;
use warnings;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use Data::Dumper;
use above 'Genome';
use Genome::Model::ClinSeq::Util qw(:all);

#Required parameters. 
my $infile = '';
my $name_column = '';
my $gene_groups = '';
my $outfile = '';

GetOptions ('infile=s'=>\$infile, 'name_column=i'=>\$name_column, 'gene_groups=s'=>\$gene_groups, 'outfile=s'=>\$outfile);

my $usage=<<INFO;

  Example usage: 
  
  geneToGeneCategories.pl  --infile=Genes.txt  --name_column=1  --gene_groups='Default'  --outfile=GenesAnnotated.txt
  
  Notes:
  This script will take an input TSV file containing a gene name column and will append columns indicating whether each gene belongs to a list of categories
  e.g. Kinases, transcription factors, etc.  These lists are pre-defined as described below.
  Whether each gene belongs to each category will be indicated with a 1 or 0.
  The order and content of lines in the input file will be maintained (including duplicate gene records)
  The order of the gene categories in your list will also be used

  Inputs:
  --infile               PATH. Any tab delimited file to be annotated containing a gene name column. **Assumed to contain a header line**.
  --name_column          INT. Number of column containing gene names (Entrez gene symbols ideally - gene name translation will be attempted).
  --gene_groups          STRING.  Name of gene group sublist (e.g. 'Default', 'Individual', 'KeggCancerPathways').  See reference annotations dir for details:
                         /gscmnt/sata132/techd/mgriffit/reference_annotations/GeneSymbolLists/config/
  --outfile              PATH.  New file with annotation columns appended.

INFO

unless ($infile && $name_column && $gene_groups && $outfile){
  print GREEN, "$usage", RESET;
  exit();
}

#Get Entrez and Ensembl data for gene name mappings
my $entrez_ensembl_data = &loadEntrezEnsemblData();

#Location of gene category files (any .txt file in this directory)
my $clinseq_annotations_dir = "/gscmnt/sata132/techd/mgriffit/reference_annotations/";
my $gene_symbol_lists_dir = $clinseq_annotations_dir . "GeneSymbolLists/";

$gene_symbol_lists_dir = &checkDir('-dir'=>$gene_symbol_lists_dir, '-clear'=>"no");
my $symbol_list_names = &importSymbolListNames('-gene_symbol_lists_dir'=>$gene_symbol_lists_dir, '-verbose'=>1);
my $master_list = $symbol_list_names->{master_list};
my @symbol_list_names = sort {$master_list->{$a}->{order} <=> $master_list->{$b}->{order}} keys %{$master_list};
my $gene_symbol_lists = &importGeneSymbolLists('-gene_symbol_lists_dir'=>$gene_symbol_lists_dir, '-symbol_list_names'=>\@symbol_list_names, '-entrez_ensembl_data'=>$entrez_ensembl_data, '-verbose'=>0);
my $master_group_list = $symbol_list_names->{master_group_list};
my $sublists = $symbol_list_names->{sublists};

unless (defined($sublists->{$gene_groups})){
  print RED, "\n\nCould not find a gene sublist with the name: $gene_groups\n\n", RESET;
  exit();
}
my $target_groups = $sublists->{$gene_groups}->{groups};

#Import gene names from input file:
my %lines;
open (IN, "$infile") || die "\n\nCould not open input file: $infile\n\n";
my $header = 1;
my $header_line = "";
my $l = 0;
while(<IN>){
  chomp($_);
  my @line = split("\t", $_);
  if ($header){
    $header = 0;
    $header_line = $_;
    next();
  }
  $l++;
  my $gene_name = $line[$name_column-1];
  $lines{$l}{record} = $_;
  $lines{$l}{gene_name} = $gene_name;

  #Fix gene names
  my $mapped_gene_name = &fixGeneName('-gene'=>$gene_name, '-entrez_ensembl_data'=>$entrez_ensembl_data, '-verbose'=>0);
  $lines{$l}{mapped_gene_name} = $mapped_gene_name;
}
close(IN);


#Intersect the gene names in the input list with each symbol list
my %member_list;
foreach my $l (keys %lines){
  my $mapped_gene_name = $lines{$l}{mapped_gene_name};
  foreach my $group_name (sort {$target_groups->{$a}->{order} <=> $target_groups->{$b}->{order}} keys %{$target_groups}){
    my @group_members = @{$master_group_list->{$group_name}->{members}};
    foreach my $member (@group_members){
      $member_list{$member}=1;
      my $symbols = $gene_symbol_lists->{$member}->{symbols};
      if ($gene_symbol_lists->{$member}->{symbols}->{$mapped_gene_name}){
        $lines{$l}{$member} = 1;
      }else{
        $lines{$l}{$member} = 0;
      }
    }
  }
}
my @member_list = keys %member_list;
my @member_list_sort = sort @member_list;


#Print output file with annotations and new column headers appended
open (OUT, ">$outfile") || die "\n\nCould not open output file for writing: $outfile\n\n";
my $append_header = join("\t", @member_list_sort);
print OUT "$header_line\tMappedGeneName\t$append_header\n";

foreach my $l (sort {$a <=> $b} keys %lines){
  my $mapped_gene_name = $lines{$l}{mapped_gene_name};
  my @append;
  foreach my $member (@member_list_sort){
    push(@append, $lines{$l}{$member});
  }
  my $append_string = join("\t", @append);
  print OUT "$lines{$l}{record}\t$mapped_gene_name\t$append_string\n";
}
close (OUT);
print "\n\n";


exit();

