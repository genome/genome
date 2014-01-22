#!/usr/bin/env genome-perl
#Written by Malachi Griffith
#For a group of ClinSeq models list basic info like build dirs

#Input:
#A list of Clinseq builds, models, or a Clinseq model-group

use strict;
use warnings;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use Data::Dumper;
use Statistics::Descriptive;
use above 'Genome';
use Genome::Model::ClinSeq::Util qw(:all);
use Genome::Model::ClinSeq::OriginalScripts::Converge qw(:all);

my $build_ids = '';
my $model_ids = '';
my $model_group_id = '';
my $amp_cutoff = '';
my $del_cutoff = '';
my $outdir = '';
my $verbose = 0;

GetOptions ('build_ids=s'=>\$build_ids, 'model_ids=s'=>\$model_ids, 'model_group_id=s'=>\$model_group_id, 
            'verbose=i'=>\$verbose);


my $usage=<<INFO;
Example usage: 

  listModelProperties.pl  --model_group_id='31779'  --verbose=1

  Specify *one* of the following as input (each model/build should be a ClinSeq model)
  --build_ids                Comma separated list of specific build IDs
  --model_ids                Comma separated list of specific model IDs
  --model_group_id           A singe genome model group ID

  --verbose                  More descriptive stdout messages

  Test Clinseq model groups:
  31779                      ClinSeq - LUCs - v1

INFO

unless (($build_ids || $model_ids || $model_group_id)){
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
if ($verbose){ print "\n\n";}

my %mb = %{$models_builds->{cases}};
foreach my $c (keys %mb){
  my $b = $mb{$c}{build};
  my $m = $mb{$c}{model};
  my $model_name = $m->name;
  my $model_id = $m->id;
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

  print "$subject_name\t$final_name\t$data_directory\t$model_name\t$model_id\t$build_id\n";
}


if ($verbose){ print "\n\n";}


exit();


