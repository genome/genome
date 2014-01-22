#!/usr/bin/env genome-perl
#Written by Jason Walker, modified by Malachi Griffith
#Get a list of patient common names from the user.  
#Use the Genome API to list information about each of these patients relating to exome or other capture data sets
  

use strict;
use warnings;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use Data::Dumper;
use above 'Genome';

#Required
my $model_group = '';

GetOptions ('model_group=s'=>\$model_group);

my $usage=<<INFO;
  Example usage: 
  
  listRNAseqBuildByModelGroup.pl  --model_group='21942'

  Attempts to find RNA-seq builds for a single model group:

INFO

unless ($model_group){
  print RED, "\n\nRequired parameter missing", RESET;
  print GREEN, "\n\n$usage", RESET;
  exit(1);
}


#genome/lib/perl/Genome/ModelGroup.pm

#Get model group object
my $mg = Genome::ModelGroup->get("id"=>$model_group);
#print Dumper $mg;

#Get display name of the model group
my $display_name = $mg->__display_name__;
#print Dumper $display_name;

#Get subjects (e.g. samples) associated with each model of the model-group
my @subjects = $mg->subjects;
#print Dumper @subjects;

#Get the members of the model-group, i.e. the models
my @models = $mg->models;
#print Dumper @models;

#Cycle through the models and get their builds
foreach my $m (@models){
  my $model_name = $m->name;
  my $model_id = $m->genome_model_id;
  my $last_complete_build_id = $m->_last_complete_build_id || "NONE_FINISHED";
  my $individual_id = $m->subject->patient->id;
  my $subject_name = $m->subject->patient->name;
  my $common_name = $m->subject->patient->common_name;
  my $sample_id = $m->subject->id;
  my $sample = Genome::Sample->get($sample_id);
  #print Dumper $m->subject->patient;
  #print Dumper $m->subject;

  my $sample_name = $sample->name || "UNDEF";
  my $extraction_type = $sample->extraction_type || "UNDEF";
  my $sample_common_name = $sample->common_name || "UNDEF";
  my $tissue_desc = $sample->tissue_desc || "UNDEF";
  my $cell_type = $sample->cell_type || "UNDEF";

  my $data_directory = "NA";
  if ($m->_last_complete_build_id){
     my $b = Genome::Model::Build->get($last_complete_build_id);
     $data_directory = $b->data_directory;
     #print Dumper $b;
  }

  #print "$common_name\t$subject_name\t$sample_name\t$model_name\t$sample_common_name\t$data_directory\n";

  #LUC9	Tumor	RNAseq	/gscmnt/gc7001/info/model_data/2880886936/build116090927/alignments/accepted_hits.bam	/gscmnt/gc7001/info/model_data/2880886936/build116090927/	/gscmnt/sata420/info/model_data/2857786885/build102671028/all_sequences.fa	GRCh37-lite-build37	
  if ($m->_last_complete_build_id){
    print "$common_name\t$sample_common_name\tRNAseq\t$data_directory"."/alignments/accepted_hits.bam\t$data_directory/\t/gscmnt/sata420/info/model_data/2857786885/build102671028/all_sequences.fa\tGRCh37-lite-build37\n";
  }


}

print "\n\n";

exit();
