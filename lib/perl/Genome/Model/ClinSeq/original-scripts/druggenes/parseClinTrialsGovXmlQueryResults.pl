#!/usr/bin/env genome-perl
#Written by Obi Griffith & Malachi Griffith

#Purpose:
#This script parses the results of queryClinTrialsGovXml.pl and attempts to reduce the result down to a more reasonable size, ordered for review
#Eliminate any line gene names in black list
#Filter down to matches with cancer terms
#Require gene to be member of an existing dgidb interaction
#Could also require drug to be in existing dgidb interaction but drug names are not present in query results and would have to be parsed from xml file

use strict;
use warnings;
use Term::ANSIColor qw(:constants);
use Getopt::Long;
use Data::Dumper;
binmode(STDOUT, ":utf8");
use Cwd 'abs_path';
use above 'Genome';
use Genome::Model::ClinSeq::Util qw(:all);
use XML::Simple;

my $lib_dir;
BEGIN{
  if (abs_path($0) =~ /(.*\/).*\/.*\.pl/){
    $lib_dir = $1."/drug-genes";
  }
}
use lib $lib_dir;
use utility qw(:all);

my $query_results_file = '';
my $kw_blacklist_file = '';
my $outfile = '';
my $xml_dir = '';
my $verbose = '';

GetOptions ('query_results_file=s'=>\$query_results_file, 'kw_blacklist_file=s'=>\$kw_blacklist_file, 'outfile=s'=>\$outfile, 'xml_dir=s'=>\$xml_dir, 'verbose=i'=>\$verbose);

my $usage=<<INFO;

  Example usage: 
  
  parseClinTrialsGovXmlQueryResults.pl  --query_results_file=/gscuser/ogriffit/Projects/DruggableGenes/KnownDruggable/ClinicalTrialsGov/interactions/query_results_summary.txt  --kw_blacklist_file=/gscmnt/sata132/techd/mgriffit/DruggableGenes/KnownDruggable/ClinicalTrialsGov/keyword_blacklist.txt  --xml_dir=/gscmnt/sata132/techd/mgriffit/DruggableGenes/KnownDruggable/ClinicalTrialsGov/xml_records_singles/  --outfile=/gscuser/ogriffit/Projects/DruggableGenes/KnownDruggable/ClinicalTrialsGov/interactions/query_results_summary_cancer_ranksorted.txt

  Details:
  --query_results_file        PATH.  File containing results of queryClinTrialsGovXml.pl
  --kw_blacklist_file         PATH.  File containing keywords that will not be allowed when attempting to match records to Gene Names
  --xml_dir                   PATH.  Directory containing all corresponding clinical trials XML files (one per record)
  --outfile                   PATH.  File to store filtered results
  --verbose                   INT.   More output

INFO

unless ($query_results_file && $kw_blacklist_file && $xml_dir && $outfile){
  print GREEN, "\n\n$usage\n\n", RESET;
  exit();
}

#Load existing drugs and genes from dgidb interactions
#Not currently using TTD because no gene symbol field in intermediate TSV file
my %drugs; my %genes;
my $drugbank_file="/gscmnt/sata132/techd/mgriffit/DruggableGenes/TSV/DrugBank_WashU_INTERACTIONS.tsv";
my $pharmgkb_file="/gscmnt/sata132/techd/mgriffit/DruggableGenes/TSV/PharmGKB_WashU_INTERACTIONS.tsv";
my $talc_file="/gscmnt/sata132/techd/mgriffit/DruggableGenes/TSV/TALC_WashU_INTERACTIONS.tsv";

open (DRUGBANK, "$drugbank_file") or die "can not open $drugbank_file\n";
my $header=<DRUGBANK>;
while (<DRUGBANK>){
  my @data=split("\t", $_);
  $genes{$data[12]}++;
  $drugs{$data[2]}++;
}
close (DRUGBANK);

open (PHARMGKB, "$pharmgkb_file") or die "can not open $pharmgkb_file\n";
$header=<PHARMGKB>;
while (<PHARMGKB>){
  my @data=split("\t", $_);
  $drugs{$data[10]}++;
  $genes{$data[21]}++;
}
close (PHARMGKB);

open (TALC, "$talc_file") or die "can not open $talc_file\n";
$header=<TALC>;
while (<TALC>){
  my @data=split("\t", $_);
  $drugs{$data[2]}++;
  $genes{$data[1]}++;
}
close (TALC);

#Import list of blacklisted keywords
my %kw_blacklist;
open(BL, "$kw_blacklist_file") || die "\n\nCould not open $kw_blacklist_file\n\n";
while(<BL>){
  chomp($_);
  $kw_blacklist{$_} = 1;
}
close(BL);

print BLUE, "\n\nParsing query results in: $query_results_file\n\n", RESET;
open (QUERYRESULTS, "$query_results_file") || die "\n\nCould not open $query_results_file\n\n";

my %filtered_lines;
my %filtered_scores;
my $l=0;
while(<QUERYRESULTS>){
  $l++;
  my $blacklist_check=0;
  my $line=$_;
  chomp($line);
  my @data=split("\t", $line);
  my $geneid=$data[0];
  my $queryterm=$data[1];
  my $nct_id=$data[2];
  my $order=$data[3];
  my $score=$data[4];
  my $title=$data[5];
  my $condition_summary=$data[6];
  my $url=$data[7];
  unless (length($queryterm)>2){next;}
  unless ($condition_summary=~/cancer|neoplasm|tumor|tumour|carcinoma|melanoma|NSCLC|Metastases|Metastasis|Metastatic|Leukemia|Lymphoma|malignant/i){
    next; #Limit to just "cancer" related drugs for now
  }
  foreach my $blacklist (keys %kw_blacklist){
    if (uc($queryterm) eq uc($blacklist)){
      $blacklist_check=1;
      last;
    }
  }
  if ($blacklist_check==1){next;} #Eliminate problematic "blacklist" gene names
  unless ($genes{uc($queryterm)}){next;} #Require gene to be in existing dgidb interaction
  

  my $drugs=&getDrugsFromXML($nct_id);


  $filtered_lines{$l}=$line;
  $filtered_scores{$l}=$score;
}
close(QUERYRESULTS);

open (OUTFILE, ">$outfile") or die "\n\ncould not write to $outfile\n\n";
foreach my $line (sort {$filtered_scores{$b} <=> $filtered_scores{$a}} keys %filtered_scores ){
  print OUTFILE "$filtered_lines{$line}\n";
}
close(OUTFILE);

sub getDrugsFromXML{
  my $rid=shift;
  my $xml_path="$xml_dir"."$rid".".xml";
  print "parsing $xml_path\n";
  my $xs1 = XML::Simple->new();

  if (-e $xml_path && -s $xml_path){
    if ($verbose){
      print CYAN, "\n\t$rid: parsing ($xml_path)\n", RESET;
    }
  }elsif(-e $xml_path && -z $xml_path){
    if ($verbose){
      print YELLOW, "\n\t$rid: empty file ($xml_path)\n", RESET;
    }
    next();
  }else{
    print RED, "\n\nCould not find XML file for record ($rid): $xml_path\n\n", RESET;
    exit();
  }

  my $xml = $xs1->XMLin($xml_path, KeyAttr => ['nct_id', 'id'] );
  my $intervention = $xml->{'intervention'};
  my @intervention_names = @{&parseXmlTree('-ref'=>$intervention, '-value_name'=>'intervention_name')};
  my $in_count = scalar(@intervention_names);
  my @intervention_types = @{&parseXmlTree('-ref'=>$intervention, '-value_name'=>'intervention_type')};
  my $it_count = scalar(@intervention_types);
  #Make sure that there is equal numbers of intervention types and names. They should correspond directly
  unless ($in_count == $it_count){
    print RED, "\n\nWarning: Interventions do not line up!\n\n", RESET;
    exit();
  }
  my $i;
  for ($i=0; $i<$it_count; $i++){
    print "$intervention_names[$i]\t$intervention_types[$i]\n";
  }

  print Dumper($intervention);
#  my @interventions = @{&parseXmlTree('-ref'=>$xml, '-value_name'=>'intervention')};
#  my @interventions = @{$intervention};
#  print Dumper(@interventions);
  print Dumper (@intervention_types);
  print Dumper (@intervention_names);

}

exit;
