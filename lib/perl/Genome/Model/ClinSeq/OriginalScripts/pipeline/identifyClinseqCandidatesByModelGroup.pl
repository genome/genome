#!/usr/bin/env genome-perl
#Written by Malachi Griffith

#Starting with a list of model-groups, summarize the status of WGS somatic, Exome somatic, RNAseq tumor, RNAseq tumor models for each

use strict;
use warnings;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use Data::Dumper;
use above 'Genome';

my $model_groups = '';

GetOptions ('model_groups=s'=>\$model_groups);

my $usage=<<INFO;
  Example usage: 
  
  identifyClinseqCandidatesByModelGroup.pl  --model_groups='18336,23405,27533'

  --model_groups     Comma separated list of model groups to search summarize

INFO

unless ($model_groups){
  print GREEN, "\n\n$usage\n\n", RESET;
  exit(1);
}
my @model_groups = split(",", $model_groups);

my %individuals;

foreach my $model_group_id (@model_groups){
  my $model_group = Genome::ModelGroup->get($model_group_id);
  my @models= $model_group->models;
  foreach my $model (@models){
    my $model_id = $model->id;
    my $subject = $model->subject;
    my $individual = $subject->patient;
    my $common_name = $individual->common_name || "NA";

    unless ($individual->isa("Genome::Individual")) {
      print "\n\nPatient of model: $model_id is not a Genome::Individual\n\n";
      exit();
    }

    #Get individual name and id
    my $individual_id = $individual->id;
    my $individual_name = $individual->name || "NULL";
    $individuals{$individual_id}{name} = $individual_name;
    
    if ($individual->common_name){
      $individuals{$individual_id}{common_name} = $common_name;
    }

    #What kind of model is this? What was the processing profile?
    my $pp_id = $model->processing_profile_id;
    my $pp = Genome::ProcessingProfile->get($pp_id);
    my $pp_name = $pp->name;
    my $pp_type = $pp->type_name;


    #Skip models that do not have a pp type of 'somatic variation' or 'rna seq'
    unless ($pp_type =~ /somatic\s+variation|rna\s+seq/){
      next();
    }

    #Get the sample common name (e.g. 
    my $scn = $subject->common_name || "NA";

    #Store models according to their types.  Use pp_name to determine whether the analysis was Exome or WGS
    my $model_type = '';
    if ($pp_type =~ /somatic\s+variation/){
      if ($pp_name =~ /wgs/i){
        $model_type = "wgs_somatic";
      }elsif($pp_name =~ /exome/i){
        $model_type = "exome_somatic";
      }else{
        print "\n\nWarning was not able to determine WGS or Exome for model $model_id from pp_name: $pp_name";
      }
    }elsif ($pp_type =~ /rna\s+seq/){
      $model_type = "rna_seq";
    }
    
    #print "\n$common_name\t$pp_id\t$pp_name\t$pp_type\t$model_type";
    if (defined($individuals{$individual_id}{$model_type})){
      my $models = $individuals{$individual_id}{$model_type};
      $models->{$model_id}->{pp_name} = $pp_name;
      $models->{$model_id}->{pp_id} = $pp_id;
      $models->{$model_id}->{model_type} = $model_type;
    }else{
      my %tmp;
      $tmp{$model_id}{pp_name} = $pp_name;
      $tmp{$model_id}{pp_id} = $pp_id;
      $tmp{$model_id}{model_type} = $model_type;
      $individuals{$individual_id}{$model_type} = \%tmp;
    }
  }
}

print Dumper %individuals;

my @model_types = qw (wgs_somatic exome_somatic rna_seq);
my %mt_counts;
foreach my $model_type (@model_types){
  $mt_counts{$model_type} = 0;
}

my $triple_threat_count = 0;
foreach my $individual_id (sort keys %individuals){

  my $common_name = "Unknown";
  if (defined($individuals{$individual_id}{common_name})){
    $common_name = $individuals{$individual_id}{common_name};
  }
  my $individual_name = $individuals{$individual_id}{name};
 
  print "\n\n$individual_id\t$individual_name\t$common_name";
  my $model_type_count = 0;
  foreach my $model_type (@model_types){
    my $model_count = 0;
    if (defined($individuals{$individual_id}{$model_type})){
      my $models = $individuals{$individual_id}{$model_type};
      $model_count = keys %{$models};
    }
    print "\n\t$model_type\t$model_count";
    if ($model_count > 0){
      $mt_counts{$model_type}++;
      $model_type_count++;
    }
  }
  if ($model_type_count >= 3){
    $triple_threat_count++;
  }
  print "\n\tModel type count = $model_type_count";
}

my $individual_count = keys %individuals;

print "\n\nFound $individual_count unique individuals (by id)";
foreach my $model_type (@model_types){
  my $count = $mt_counts{$model_type};
  print "\n\t$model_type = $count";
}

print "\n\nIndividuals with all three data types = $triple_threat_count\n\n";

exit();






