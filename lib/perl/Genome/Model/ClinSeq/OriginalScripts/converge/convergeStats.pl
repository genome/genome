#!/usr/bin/env genome-perl
#Written by Malachi Griffith
#For a group of ClinSeq models, converge all Stats.tsv files (RNA-seq library quality metrics, SNV concordance, etc.)
#All stats should be pre-calculated
#Produce an output matrix in which each row is metric and each column is a library

#Input:
#A list of Clinseq builds, models, or a Clinseq model-group

use strict;
use warnings;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use Data::Dumper;
use above 'Genome';
use Genome::Model::ClinSeq::Util qw(:all);
use Genome::Model::ClinSeq::OriginalScripts::Converge qw(:all);

my $build_ids = '';
my $model_ids = '';
my $model_group_id = '';
my $outfile = '';
my $verbose = 0;

GetOptions ('build_ids=s'=>\$build_ids, 'model_ids=s'=>\$model_ids, 'model_group_id=s'=>\$model_group_id, 
            'outfile=s'=>\$outfile, 'verbose=i'=>\$verbose);

my $usage=<<INFO;
  Example usage: 

  convergeStats.pl  --model_group_id='32264'  --outfile=TechD_RNAseq_Metrics.tsv  --verbose=1

  Specify *one* of the following as input (each model/build should be a ClinSeq model)
  --build_ids                Comma separated list of specific build IDs
  --model_ids                Comma separated list of specific model IDs
  --model_group_id           A singe genome model group ID

  --outfile                  Path of the output file to be written
  --verbose                  More descriptive stdout messages

  Test Clinseq model groups:
  32264                      ClinSeq - TechD RNA-seq library comparison - SPIA trimmed data
  32181                      ClinSeq - TechD RNA-seq library comparison
  30176                      ClinSeq - RNA-seq (polyA) vs cDNA Capture (NimbleGen Exome v2) Evaluation

INFO

unless (($build_ids || $model_ids || $model_group_id) && $outfile){
  print RED, "\n\nRequired parameter missing", RESET;
  print GREEN, "\n\n$usage", RESET;
  exit(1);
}

#Get the models/builds
if ($verbose){print BLUE, "\n\nGet genome models/builds for supplied list", RESET;}
my $models_builds;
if ($build_ids){
  my @build_ids = split(",", $build_ids);
  $models_builds = &getModelsBuilds('-builds'=>\@build_ids, '-verbose'=>$verbose);
}elsif($model_ids){
  my @model_ids = split(",", $model_ids);
  $models_builds = &getModelsBuilds('-models'=>\@model_ids, '-verbose'=>$verbose);
}elsif($model_group_id){
  $models_builds = &getModelsBuilds('-model_group_id'=>$model_group_id, '-verbose'=>$verbose);
}else{
  print RED, "\n\nCould not obtains models/builds - check input to convergeCufflinksExpression.pl\n\n", RESET;
  exit();
}

my $target_file_name = "Stats.tsv";
my %files = %{&getStatsFiles('-models_builds'=>$models_builds, '-verbose'=>$verbose)};

#Build a hash of possible stats values.  Key it on:
#Question + Data_Type + Analysis_Type + Statistic_Type + Source_File - Make sure these keys are unique within a patient (column_name)!!
my $result = &parseMetrics('-files'=>\%files);
my $metrics = $result->{'metrics'};
my $metric_list = $result->{'metric_list'};

#Build the header line
my @column_names;
foreach my $build_id (sort {$files{$a}->{column_name} cmp $files{$b}->{column_name}} keys %files){
  my $column_name = $files{$build_id}{column_name};
  push(@column_names, $column_name);
}
my $column_names_s = join("\t", @column_names);
my $header_line = "Question\tData_Type\tAnalysis_Type\tStatistic_Type\tExtra_Description\tFile_Source\t$column_names_s";

#No go through each build and print out the summary
open (OUT, ">$outfile") || die "\n\nCould not open output file: $outfile\n\n";
print OUT "$header_line\n";
foreach my $metric (sort {$metric_list->{$a}->{order} <=> $metric_list->{$b}->{order}} keys %{$metric_list}){
  my $metric_info = $metric_list->{$metric}->{question} . "\t" . $metric_list->{$metric}->{data_type} . "\t" . $metric_list->{$metric}->{analysis_type} . "\t" . $metric_list->{$metric}->{statistic_type} . "\t" . $metric_list->{$metric}->{extra_description} . "\t" . $metric_list->{$metric}->{distinct_name};

  my @build_metrics;

  foreach my $build_id (sort {$files{$a}->{column_name} cmp $files{$b}->{column_name}} keys %files){
    #Watch for undefined metrics
    my $answer = "NA";
    if (defined($metrics->{$build_id}->{$metric})){
      $answer = $metrics->{$build_id}->{$metric}->{answer};
    }
    push(@build_metrics, $answer);
  }
  my $build_metrics_s = join("\t", @build_metrics);

  print OUT "$metric_info\t$build_metrics_s\n";
}
close(OUT);
print "\n\nWrote output to: $outfile\n\n";

exit();




#######################################################################################################################
#Get the desired files from each model                                                                                #
#######################################################################################################################
sub getStatsFiles{
  my %args = @_;
  my $models_builds = $args{'-models_builds'};
  my $verbose = $args{'-verbose'};

  my %files;
  my $fc = 0;
  if ($verbose){print BLUE, "\n\nGet all Stats files within these builds that match 'Stats.tsv'", RESET;}
  my %mb = %{$models_builds->{cases}};
  foreach my $c (keys %mb){
    my $b = $mb{$c}{build};
    my $m = $mb{$c}{model};
    my $model_name = $m->name;
    my $data_directory = $b->data_directory;
    my $subject_name = $b->subject->name;
    my $subject_common_name = $b->subject->common_name;
    my $build_id = $b->id;

    #If the subject name is not defined, die
    unless ($subject_name){
      print RED, "\n\nCould not determine subject name for build: $build_id\n\n", RESET;
      exit(1);
    }

    my $final_name = "Unknown";
    if ($subject_name){$final_name = $subject_name;}
    if ($subject_common_name){$final_name = $subject_common_name;}

    if ($verbose){print BLUE, "\n\t$final_name\t$build_id\t$data_directory", RESET;}

    my @files;

    #Find all 'Stats.tsv' files in the build dir
    my $find_cmd = "find $data_directory -name Stats.tsv";
    if ($verbose){print YELLOW, "\n\t\t$find_cmd", RESET;}
    my @tmp = `$find_cmd`;
    chomp(@tmp);
    my $file_count = scalar(@tmp);
    my %file_info;
    foreach my $path (@tmp){
      if ($path =~ /^$data_directory\/.*?\/(.*)/){
        my $distinct_name = $1;
        $file_info{$path}{distinct_name} = $distinct_name;
      }else{
        print RED, "\n\nCould not determine distinct name from path\n$path\n\n", RESET;
        exit(1);
      }
    }

    if ($file_count > 0){
      $files{$build_id}{stats_files} = \%file_info;
      $files{$build_id}{final_name} = $final_name;
      $files{$build_id}{subject_name} = $subject_name;
      $files{$build_id}{model_name} = $model_name;
    }else{
      if ($verbose){print YELLOW, "\n\tWarning.  Could not find any Stats.tsv files for build: $build_id", RESET;}
    }
  }

  #Determine a suitable name for data columns in the output file.  It must be distinct accross the builds.  Try these combinations in order:
  # final_name
  # final_name + subject_name
  # model_name
  # final_name + subject_name + build_id
  my $file_count =  keys %files;
  unless ($file_count > 1){
    print RED, "\n\nFound $file_count files to join (need at least 2 ...)\n\nLooked for files of name: Stats.tsv\n\n", RESET;
    exit(1);
  }
  my %finalnames;
  my %finalnames_subjectnames;
  my %modelnames;
  my %finalnames_subjectnames_builds;
  foreach my $bid (keys %files){
    my $final_name = $files{$bid}{final_name};
    my $subject_name = $files{$bid}{subject_name};
    my $model_name = $files{$bid}{model_name};
    my $id1 = $final_name;
    my $id2 = "$final_name"."_"."$subject_name";
    my $id3 = $model_name;
    my $id4 = "$final_name"."_"."$subject_name"."_"."$bid";
    
    $finalnames{$id1}=1;
    $finalnames_subjectnames{$id2}=1;
    $modelnames{$id3}=1;
    $finalnames_subjectnames_builds{$id4}=1;

    $files{$bid}{id1} = $id1;
    $files{$bid}{id2} = $id2;
    $files{$bid}{id3} = $id3;
    $files{$bid}{id4} = $id4;
  }

  my $id_choice;
  if ($file_count == keys %finalnames){
    $id_choice = "id1";
  }elsif($file_count == keys %finalnames_subjectnames){
    $id_choice = "id2";
  }elsif($file_count == keys %modelnames){
    $id_choice = "id3";
  }elsif($file_count == keys %finalnames_subjectnames_builds){
    $id_choice = "id4";
  }else{
    print RED, "\n\nCould not generate a suitable set of distinct labels for the data files to be joined...\n\n", RESET;
    exit();
  }

  foreach my $fc (keys %files){
    $files{$fc}{column_name} = $files{$fc}{"$id_choice"};
  }

  return(\%files);
}


#######################################################################################################################
#Build a hash of metrics across the input builds and their stats files                                                #
#######################################################################################################################
sub parseMetrics{
  my %args = @_;
  my %files = %{$args{'-files'}};

  my %grand_metrics;
  my %metric_list;
  my $o = 0;
  foreach my $build_id (sort keys %files){
    my %local_metrics;
    my %stats_files = %{$files{$build_id}{stats_files}};
    my $column_name = $files{$build_id}{column_name};
    foreach my $file (sort keys %stats_files){
      my $distinct_name = $stats_files{$file}{distinct_name};

      open (STATS, "$file") || die "\n\nCould not open stats file: $file\n\n";
      my %columns;
      my $header = 1;
      while(<STATS>){
        chomp($_);
        my @line = split("\t", $_);
        
        #Perform sanity check on header to make sure it conforms to the ClinSeq stats.tsv standard
        if ($header){
          my $p = 0;
          foreach my $col (@line){
            $columns{$col}{position} = $p;
            $p++;
          }
          $header = 0;
          unless ($columns{'Question'} && $columns{'Answer'} && $columns{'Data_Type'} && $columns{'Analysis_Type'} && $columns{'Statistic_Type'} && $columns{'Extra_Description'}){
            print RED, "\n\nRequired column missing from file:\nRequired columns: Question, Answer, Data_Type, Analysis_Type, Statistic_Type, Extra_Description\nFile: $file\nHeader: @line\n\n", RESET;
            exit(1);
          }
          next();
        }

        #Parse the metrics data and store in a hash keyed on: build id (one per column in the final output) & a concatenated string unique to the metric
        my $question = $line[$columns{'Question'}{position}];
        my $answer = $line[$columns{'Answer'}{position}];
        my $data_type = $line[$columns{'Data_Type'}{position}];
        my $analysis_type = $line[$columns{'Analysis_Type'}{position}];
        my $statistic_type = $line[$columns{'Statistic_Type'}{position}];
        my $extra_description = $line[$columns{'Extra_Description'}{position}];

        my $metric_key = $question.$data_type.$analysis_type.$statistic_type.$extra_description.$distinct_name;

        #Clean up the key a bit
        $metric_key =~ s/\'//g;
        $metric_key =~ s/\%//g;
        $metric_key =~ s/\=//g;

        #A metric must not be duplicated within a metrics file!
        if (defined($local_metrics{$metric_key})){
          print RED, "\n\nMetric on this line appears to be a duplicate:\n@line\n\n$file\n\n", RESET;
          exit(1);
        }else{
          $local_metrics{$metric_key} = 1;
          
          #Store the descriptive metric info once for all builds in a distinct list
          unless (defined($metric_list{$metric_key})){
            $o++;
            $metric_list{$metric_key}{order} = $o;
            $metric_list{$metric_key}{question} = $question;
            $metric_list{$metric_key}{data_type} = $data_type;
            $metric_list{$metric_key}{analysis_type} = $analysis_type;
            $metric_list{$metric_key}{statistic_type} = $statistic_type;
            $metric_list{$metric_key}{extra_description} = $extra_description;
            $metric_list{$metric_key}{distinct_name} = $distinct_name;
          }
          #Store the actual metric values for each individual build
          $grand_metrics{$build_id}{$metric_key}{answer} = $answer;

        }
      }
      close(STATS);
    }
  }

  my %results;
  $results{metrics} = \%grand_metrics;
  $results{metric_list} = \%metric_list;

  return(\%results);
}

