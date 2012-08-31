#!/usr/bin/env genome-perl
#Written by Malachi Griffith
#Starting with a model-group for any combination of the following: WGS somatic variation model, Exome somatic variation model, RNA-seq model
#Identify those models that can be grouped by individual and create clinseq models for these:
#The rna-seq model group may contain both 'tumor' and 'normal'

use strict;
use warnings;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use Data::Dumper;
use above 'Genome';

my $wgs_somatic_model_group = '';
my $exome_somatic_model_group = '';
my $rnaseq_model_group = '';
my $normal_tissue_names = '';
my $project_label = '';
my $outfile = '';


GetOptions ('wgs_somatic_model_group=s'=>\$wgs_somatic_model_group, 'exome_somatic_model_group=s'=>\$exome_somatic_model_group, 'rnaseq_model_group=s'=>\$rnaseq_model_group, 'normal_tissue_names=s'=>\$normal_tissue_names, 'project_label=s'=>\$project_label, 'outfile=s'=>\$outfile);

my $usage=<<INFO;
  Example usage: 
  
  createClinseqJobs.pl  --wgs_somatic_model_group=25481  --rnaseq_model_group=29586  --outfile=ALS_Clinseq_Jobs.txt

  --wgs_somatic_model_group     Group of WGS somatic variation models (name or ID)
  --exome_somatic_model_group   Group of Exome somatic variation models (name or ID)
  --rnaseq_model_group          Group of RNAseq models (name or ID)
  --normal_tissue_names         Define a list of names that could be used to specify 'normal' tissue.  All others will be considered tumor.
                                e.g. 'blood|skin', 'occipital cortex|blood'
  --project_label               A label for the project that will be included in the clinseq model names (e.g. ALS, spinal cord vs. occiptal cortex)
  --outfile                     Output file where commands will be stored

INFO

unless (($wgs_somatic_model_group || $exome_somatic_model_group || $rnaseq_model_group) && $outfile){
  print GREEN, "\n\n$usage\n\n", RESET;
  exit(1);
}

my %individuals;
my @ntn = split("\\|", $normal_tissue_names);
my %normal_tissue_names;
foreach my $ntn (@ntn){
  $normal_tissue_names{$ntn}=1;
}
if (keys %normal_tissue_names > 0){
  print YELLOW, "\n\nThe normal tissue names ($normal_tissue_names) will be used to identify 'normal' RNAseq models, the rest will be considered tumors\n\n", RESET;
}


if ($wgs_somatic_model_group){
  my $mg = Genome::ModelGroup->get($wgs_somatic_model_group);
  my @models= $mg->models;
  foreach my $model (@models){
    my $model_id = $model->id;
    my $sample = $model->subject;
    my $sample_name = $sample->name;
    my $individual = $sample->patient;
    my $individual_name = $individual->name;
    my $tissue_desc = $sample->tissue_desc || "NA";
    my $common_name = $individual->common_name || $sample->name;
    #print "\nWGS\t$model_id\t$individual_name\t$common_name\t$tissue_desc";
    $individuals{$individual_name}{wgs_somatic_model_id} = $model_id;
    $individuals{$individual_name}{common_name} = $common_name;
  }
}
if ($exome_somatic_model_group){
  my $mg = Genome::ModelGroup->get($exome_somatic_model_group);
  my @models= $mg->models;
  foreach my $model (@models){
    my $model_id = $model->id;
    my $sample = $model->subject;
    my $sample_name = $sample->name;
    my $individual = $sample->patient;
    my $individual_name = $individual->name;
    my $tissue_desc = $sample->tissue_desc || "NA";
    my $common_name = $individual->common_name || $sample->name;
    #print "\nExome\t$model_id\t$individual_name\t$common_name\t$tissue_desc";
    $individuals{$individual_name}{exome_somatic_model_id} = $model_id;
    $individuals{$individual_name}{common_name} = $common_name;
  }
}
if ($rnaseq_model_group){
  my $mg = Genome::ModelGroup->get($rnaseq_model_group);
  my @models= $mg->models;
  foreach my $model (@models){
    my $model_id = $model->id;
    my $sample = $model->subject;
    my $sample_name = $sample->name;
    my $individual = $sample->patient;
    my $individual_name = $individual->name;
    my $tissue_desc = $sample->tissue_desc || "NA";
    my $common_name = $individual->common_name || $sample->name;
    #print "\nRNAseq\t$model_id\t$individual_name\t$common_name\t$tissue_desc";

    #Set tumor and normal RNAseq models if available - use '$normal_tissue_desc' to define which is 'normal'.  All others will be considered 'tumor' (e.g. 'occipital cortex')
    if (defined($normal_tissue_names{$tissue_desc})){
      $individuals{$individual_name}{normal_rnaseq_model_id} = $model_id;
      $individuals{$individual_name}{common_name} = $common_name;
    }else{
      $individuals{$individual_name}{tumor_rnaseq_model_id} = $model_id;
      $individuals{$individual_name}{common_name} = $common_name;
    }
  }
}

#genome model define clin-seq  --model-name='ClinSeq - ALSXXX - Spinal vs Occipital/Blood - (Nov. 2011 PP) '  --processing-profile='November 2011 Clinical Sequencing'  --wgs-model=''  --tumor-rnaseq-model=''  --normal-rnaseq-model=''
open (OUT, ">$outfile") || die "\n\nCould not open output file: $outfile\n\n";
foreach my $individual (sort keys %individuals){
  my ($wgs_somatic_model_id, $exome_somatic_model_id, $tumor_rnaseq_model_id, $normal_rnaseq_model_id, $common_name) = 0;
  $common_name = $individuals{$individual}{common_name};

  my $model_name = "ClinSeq - $individual ($common_name)";
  if ($project_label){
    $model_name .= " - $project_label";
  }

  my $cmd = "#$individual ($common_name)\ngenome model define clin-seq  --model-name='$model_name'  --processing-profile='November 2011 Clinical Sequencing'";
  if (defined($individuals{$individual}{wgs_somatic_model_id})){
    $cmd .= "  --wgs-model='$individuals{$individual}{wgs_somatic_model_id}'";
  }
  if (defined($individuals{$individual}{exome_somatic_model_id})){
    $cmd .= "  --exome-model='$individuals{$individual}{exome_somatic_model_id}'";
  }
  if (defined($individuals{$individual}{tumor_rnaseq_model_id})){
    $cmd .= "  --tumor-rnaseq-model='$individuals{$individual}{tumor_rnaseq_model_id}'";
  }
  if (defined($individuals{$individual}{normal_rnaseq_model_id})){
    $cmd .= "  --normal-rnaseq-model='$individuals{$individual}{normal_rnaseq_model_id}'";
  }
  $cmd .= "\ngenome model build start ''\n\n";
  print OUT "$cmd";
}
close(OUT);

print BLUE, "\nWrote commands to: $outfile\n\n", RESET;


exit();






