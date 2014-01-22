#!/usr/bin/env genome-perl
use strict;
use warnings;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use Data::Dumper;

#Import the SNV summary stats 'Stats.tsv' file for all LUCs - produce a table compiling the stats for all patients
#Create a new list of SNVs based on the 'VariantExpressionSummary.tsv' file for each patient
#- To this file add the hg18 coordinates back on so that both hg18 and hg19 coords are provide side by side
#- Create a new file merging the VariantExpressionSummary.tsv onto the original SNV file imported for each patient

my $working_dir = "/gscmnt/sata132/techd/mgriffit/luc/rnaseq_vs_snv/";
my $readcounts_files_merge_dir_hg19 = "/gscmnt/sata132/techd/mgriffit/luc/rnaseq_vs_snv/readcounts_merge_hg19/";
opendir(DIR, "$readcounts_files_merge_dir_hg19") || die $!;
my %files;
while (my $file = readdir(DIR)) {
  chomp($file);
  if ($file =~ /LUC/){
    my $patient = $file;
    my $stats_file = $readcounts_files_merge_dir_hg19 . $patient . "/summary/Stats.tsv";
    my $processed_data_file = $readcounts_files_merge_dir_hg19 . $patient . "/summary/VariantExpressionSummary.tsv";
    $files{$patient}{stats_file} = $stats_file;
    $files{$patient}{processed_data_file} = $processed_data_file;
  }
}
closedir(DIR);

my %stats;
my %patients;
foreach my $patient (sort keys %files){
  my $po = 0;
  if ($patient =~ /LUC(\d+)/){
    $po = $1;
  }
  $patients{$patient}{order}=$po;
  my $stats_file = $files{$patient}{stats_file};
  open(STATS, "$stats_file") || die $!;
  my $header = 1;
  my $o = 0;
  while(<STATS>){
    chomp($_);
    my @line = split("\t", $_);
    if ($header){
      $header = 0;
      next();
    }
    $o++;
    $stats{$line[0]}{$patient}{order} = $o;
    $stats{$line[0]}{$patient}{answer} = $line[1];
    $stats{$line[0]}{$patient}{data_type} = $line[2];
    $stats{$line[0]}{$patient}{analysis_type} = $line[3];
    $stats{$line[0]}{$patient}{statistic_type} = $line[4];
    $stats{$line[0]}{$patient}{description} = $line[5];
  }
  close(STATS);
}

#Build a header
my @patient_list;
foreach my $patient (sort {$patients{$a}{order} <=> $patients{$b}{order}} keys %patients){
  push(@patient_list, $patient);
}
my $patient_list_string = join("\t", @patient_list);
my $header_line = "Query\tQueryOrder\t$patient_list_string\tQueryDescription\n";
print $header_line;
foreach my $question (keys %stats){
  my @values;
  my $order = 0;
  my $description = "NA";
  foreach my $patient (sort {$patients{$a}{order} <=> $patients{$b}{order}} keys %patients){
    push(@values, $stats{$question}{$patient}{answer});
    $order = $stats{$question}{$patient}{order};
    $description = $stats{$question}{$patient}{description};
  }
  my $values_string = join("\t", @values);
  my $data_line = "$question\t$order\t$values_string\t$description\n";
  print $data_line;
}




exit();





