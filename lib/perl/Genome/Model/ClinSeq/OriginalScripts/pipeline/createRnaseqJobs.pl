#!/usr/bin/env genome-perl
#Written by Malachi Griffith
#From a list of patient common names, individual names, sample_names, or even a model group:
#1.) Get a list of RNA samples for those
#2.) Summarize what was found
#3.) Create commands to generate and launch RNAseq models for each unique sample

use strict;
use warnings;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use Data::Dumper;
use above 'Genome';

#Required
my $common_names = '';
my $individual_names = '';
my $sample_names = '';
my $model_groups = '';
my $processing_profile = '';
my $reference_sequence_build_id = '';
my $annotation_build = '';
my $outfile = '';

GetOptions ('common_names=s'=>\$common_names, 'individual_names=s'=>\$individual_names, 'sample_names=s'=>\$sample_names, 'model_groups=s'=>\$model_groups, 
            'processing_profile=s'=>\$processing_profile, 'reference_sequence_build_id=s'=>\$reference_sequence_build_id, 'annotation_build=s'=>\$annotation_build,
            'outfile=s'=>\$outfile);

my $usage=<<INFO;
  Example usage: 
  
  createRnaseqJobs.pl  --processing_profile='2754795'  --reference_sequence_build_id='106942997'  --annotation_build='124434505'  --outfile=LUC_RNAseq_Jobs.txt  --common_names='LUC1,LUC2'
  createRnaseqJobs.pl  --processing_profile='2754795'  --reference_sequence_build_id='106942997'  --annotation_build='124434505'  --outfile=LUC_RNAseq_Jobs.txt  --individual_names='H_JG-379,H_JG-990'
  createRnaseqJobs.pl  --processing_profile='2754795'  --reference_sequence_build_id='106942997'  --annotation_build='124434505'  --outfile=LUC_RNAseq_Jobs.txt  --sample_names='H_JG-990-S.2514,H_JG-379-1029757'  
  createRnaseqJobs.pl  --processing_profile='2754795'  --reference_sequence_build_id='106942997'  --annotation_build='124434505'  --outfile=LUC_RNAseq_Jobs.txt  --model_groups='21942'

  Create a series of RNAseq model creation and building jobs starting with a list of common names, individual names, or sample names
  Note that you can specify any combination of the following and a distinct list of samples to work on will be obtained
  --common_names                  Comma separated list of patient common names
  --individual_names              Comma separated list of individual names
  --sample_names                  Comma separated list of sample names
  --model_groups                  Comma separated list of model-groups for models involving RNA-samples
  --processing_profile            ID or full name of the desired processing profile
  --reference_sequence_build_id   ID of the reference sequence build to use
  --annotation_build              ID of the reference annotation build to use

  --outfile                       Name of file where commands and notes will be written


INFO

unless (($common_names || $individual_names || $sample_names || $model_groups) && $processing_profile && $reference_sequence_build_id && $annotation_build && $outfile){
  print GREEN, "\n\n$usage\n\n", RESET;
  exit(1);
}

#Get the processing profile ID/name - starting with a processing_profile name or ID
my $pp = Genome::ProcessingProfile->get($processing_profile);
my $pp_id;
my $pp_name;
unless($pp){
  $pp = Genome::ProcessingProfile->get(name=>$processing_profile);
}
if ($pp){
  $pp_id = $pp->id;
  $pp_name = $pp->name;
  print BLUE, "\n\nUsing processing profile: $pp_id ($pp_name)", RESET;
}else{
  print RED, "\n\nCould not get processing-profile object from user specified: $processing_profile\n\n", RESET;
  exit(1);
}

#Get the reference sequence build info and remind the user what you are going to use
my $reference_model_build = Genome::Model::Build->get($reference_sequence_build_id);
my $reference_model_name = $reference_model_build->name;
my $reference_model_display_name = $reference_model_build->__display_name__;
if ($reference_model_build && $reference_model_display_name){
  print BLUE, "\n\nUsing reference annotation build: $reference_model_display_name", RESET;
}else{
  print RED, "\n\nCould not get reference sequence build id from user specified: $reference_sequence_build_id\n\n", RESET;
  exit(1);
}

#Check input lists
print BLUE, "\n\nParsing input list(s)", RESET;
my @common_names;
my @individual_names;
my @sample_names;
my @model_groups;

my $grand_count = 0;
if ($common_names){
  @common_names = split(",", $common_names);
  my $count = scalar(@common_names);
  $grand_count += $count;
  if ($count > 0){
    print BLUE, "\n\tFound $count common_names: @common_names", RESET;
  }
}
if ($individual_names){
  @individual_names = split(",", $individual_names);
  my $count = scalar(@individual_names);
  $grand_count += $count;
  if ($count > 0){
    print BLUE, "\n\tFound $count individual_names: @individual_names", RESET;
  }
}
if ($sample_names){
  @sample_names = split(",", $sample_names);
  my $count = scalar(@sample_names);
  $grand_count += $count;
  if ($count > 0){
    print BLUE, "\n\tFound $count sample_names: @sample_names", RESET;
  }
}
if ($model_groups){
  @model_groups = split(",", $model_groups);
  my $count = scalar(@model_groups);
  $grand_count += $count;
  if ($count > 0){
    print BLUE, "\n\tFound $count model_groups: @model_groups", RESET;
  }
}
print BLUE, "\nFound a grand total of $grand_count input items", RESET;


#Get a list of sample objects from the starting list of common_names, individuals, or samples
print BLUE, "\n\nSearching for samples associated with input list(s)", RESET;
my @aggregate_samples;
if ($common_names){
  print BLUE, "\n\tBy common name(s)", RESET;
  foreach my $cn (@common_names){
    my $individual = Genome::Individual->get(common_name => $cn);
    unless ($individual){
      print YELLOW, "\n\t\tCould not obtain Genome::Individual object from: $cn", RESET;
      next();
    }
    my $individual_name = $individual->name;
    my @samples = $individual->samples;
    foreach my $sample (@samples) {
      my $extraction_type = $sample->extraction_type;
      my $sample_name = $sample->name;
      my $individual = $sample->patient;
      my $common_name = $individual->common_name || "unknown";
      my $tissue_desc = $sample->tissue_desc || "unknown";
      my $cell_type = $sample->cell_type || "unknown";
      my $sample_common_name = $sample->common_name || "unknown";
      if ($extraction_type){
        if ($extraction_type eq "rna"){
          print BLUE, "\n\t\t$common_name. Found RNA sample ($sample_name) for $individual_name of tissue description: $tissue_desc and cell type: $cell_type and sample common name: $sample_common_name", RESET;
          push (@aggregate_samples, $sample);
        }
      }else{
        print YELLOW, "\n\t\t$common_name. Could not find extraction_type for sample: $sample_name", RESET;
      }
    }
  }
}
if ($individual_names){
  print BLUE, "\n\tBy individual name(s)", RESET;
  foreach my $in (@individual_names){
    my $individual = Genome::Individual->get(name => $in);
    unless($individual){
      $individual = Genome::Individual->get(id => $in);
    }
    unless ($individual){
      print YELLOW, "\n\t\tCould not obtain Genome::Individual object from: $in", RESET;
      next();
    }
    my $individual_name = $individual->name;
    my @samples = $individual->samples;
    foreach my $sample (@samples) {
      my $extraction_type = $sample->extraction_type;
      my $sample_name = $sample->name;
      my $individual = $sample->patient;
      my $common_name = $individual->common_name || "unknown";
      my $tissue_desc = $sample->tissue_desc || "unknown";
      my $cell_type = $sample->cell_type || "unknown";
      my $sample_common_name = $sample->common_name || "unknown";
      if ($extraction_type){
        if ($extraction_type eq "rna"){
          print BLUE, "\n\t\t$common_name. Found RNA sample ($sample_name) for $individual_name of tissue description: $tissue_desc and cell type: $cell_type and sample common name: $sample_common_name", RESET;
          push (@aggregate_samples, $sample);
        }
      }else{
        print YELLOW, "\n\t\t$common_name. Could not find extraction_type for sample: $sample_name", RESET;
      }
    }
  }
}
if ($sample_names){
  print BLUE, "\n\tBy sample name(s)", RESET;
  foreach my $sn (@sample_names){
    my $sample = Genome::Sample->get(name=>$sn);
    unless ($sample){
      $sample = Genome::Sample->get(id=>$sn);
    }
    unless ($sample){
      print YELLOW, "\n\t\tCould not obtain Genome::Sample object from: $sn", RESET;
      next();
    }
    my $individual = $sample->patient;
    my $individual_name = $individual->name;
    my $extraction_type = $sample->extraction_type;
    my $sample_name = $sample->name;
    my $common_name = $individual->common_name || "unknown";
    my $tissue_desc = $sample->tissue_desc || "unknown";
    my $cell_type = $sample->cell_type || "unknown";
    my $sample_common_name = $sample->common_name || "unknown";
    if ($extraction_type){
      if ($extraction_type eq "rna"){
        print BLUE, "\n\t\t$common_name. Found RNA sample ($sample_name) for $individual_name of tissue description: $tissue_desc and cell type: $cell_type and sample common name: $sample_common_name", RESET;
        push (@aggregate_samples, $sample);
      }
    }else{
      print YELLOW, "\n\t\t$common_name. Could not find extraction_type for sample: $sample_name", RESET;
    }
  }
}
if ($model_groups){
  print BLUE, "\n\tBy model group(s)", RESET;
  foreach my $mg (@model_groups){
    my $model_group = Genome::ModelGroup->get($mg);
    unless($model_group) {
      $model_group = Genome::ModelGroup->get(name=>$mg);
    }
    unless($model_group) {
      print YELLOW, "\n\t\tCould not obtain Genome::ModelGroup object from: $mg", RESET;
      next();
    }
    my @models = $model_group->models;
    my @samples = sort {$a->name cmp $b->name } map {$_->subject} @models;
    foreach my $sample (@samples) {
      my $extraction_type = $sample->extraction_type;
      my $sample_name = $sample->name;
      my $individual = $sample->patient;
      my $individual_name = $individual->name;
      my $common_name = $individual->common_name || "unknown";
      my $tissue_desc = $sample->tissue_desc || "unknown";
      my $cell_type = $sample->cell_type || "unknown";
      my $sample_common_name = $sample->common_name || "unknown";
      if ($extraction_type){
        if ($extraction_type eq "rna"){
          print BLUE, "\n\t\t$common_name. Found RNA sample ($sample_name) for $individual_name of tissue description: $tissue_desc and cell type: $cell_type and sample common name: $sample_common_name", RESET;
          push (@aggregate_samples, $sample);
        }
      }else{
        print YELLOW, "\n\t\t$common_name. Could not find extraction_type for sample: $sample_name", RESET;
      }
    }
  }
}


#Go through the aggregate list of sample objects and create a list without duplicates
my @final_samples;
my %sample_list;
foreach my $sample (@aggregate_samples){
  my $id = $sample->id;
  unless ($sample_list{$id}){
    push(@final_samples, $sample);
    $sample_list{$id}=1;
  }
}
my $final_sample_count = scalar(@final_samples);
print BLUE, "\n\nFound a grand total of $final_sample_count samples to be processed", RESET;

#Generate output commands
print BLUE, "\n\nPrinting output commands and notes to: $outfile", RESET;
my %output;
my $s = 0;
foreach my $sample (@final_samples){
  $s++;
  my $extraction_type = $sample->extraction_type;
  my $sample_name = $sample->name;
  my $individual = $sample->patient;
  my $individual_name = $individual->name;
  my $common_name = $individual->common_name || "unknown";
  my $tissue_desc = $sample->tissue_desc || "unknown";
  my $cell_type = $sample->cell_type || "unknown";
  my $sample_common_name = $sample->common_name || "unknown";

  #Add a sort value to the sample objects
  my $sort_name;
  if ($individual->common_name){
    my $cn = $common_name;
    if ($cn =~ /^(\w{3})(\d{1})$/){
      $cn = "$1"."000$2";
    }elsif($cn =~ /^(\w{3})(\d{2})$/){
      $cn = "$1"."00$2";
    }elsif($cn =~ /^(\w{3})(\d{3})$/){
      $cn = "$1"."0$2";
    };
    $sort_name = $cn . "_" . $sample_common_name;
  }else{
    $sort_name = $individual_name;
  }

  my $sample_label_name = $individual_name;
  if ($individual->common_name){
    $sample_label_name = "$common_name ($individual_name)";
  }
  my $tissue_label_name = "unknown tissue";
  if ($sample->common_name){
    $tissue_label_name = $sample->common_name;
  }
  if ($sample->common_name && $sample->tissue_desc){
    $tissue_label_name = $sample->common_name . " (". $sample->tissue_desc . ")";
  }
  my $model_name = "$sample_label_name - $tissue_label_name - RNAseq - PP$pp_id - $reference_model_name";

  #Remove characters that would prevent proper parsing of the model name downstream (e.g. ',')
  $model_name =~ s/\,//g;

  #Example commands to be created
  #genome model define rna-seq  --processing-profile='Human GRCh37/Ensembl 58_37c_v2 - TopHat v1.4.0 - Cufflinks v1.3.0 - Picard v1.52 - Mask rRNA_MT - FAR Trim SPIA Allow Ns Min 25bp - reference only'  --reference-sequence-build=106942997  --subject='H_JG-3631-S.10804'  --model-name='LUC4 - PP2694793 - cDNA capture'
  #genome model instrument-data assign --model=''  --all
  #genome model build start ''


  $output{$s}{sort_name} = $sort_name;
  $output{$s}{header} = "\n#Common_name = $common_name, Individual_name = $individual_name, Sample_name = $sample_name, Sample_common_name = $sample_common_name, Cell_type = $cell_type, Tissue_desc = $tissue_desc";
  $output{$s}{clinseq_cmd} = "\ngenome model define rna-seq  --processing-profile='$pp_id'  --reference-sequence-build='$reference_sequence_build_id'  --annotation-build='$annotation_build'  --subject='$sample_name'  --model-name='$model_name'";
  $output{$s}{id_cmd} = "\ngenome model instrument-data assign --model='$model_name'  --all";
  $output{$s}{build_cmd} = "\ngenome model build start '$model_name'\n";
}

open (OUT, ">$outfile") || die "\n\nCould not open output file: $outfile\n\n";
foreach my $s (sort {$output{$a}{sort_name} cmp $output{$b}{sort_name}} keys %output){
  print OUT $output{$s}{header}.$output{$s}{clinseq_cmd}.$output{$s}{id_cmd}.$output{$s}{build_cmd};
}
close(OUT);


print "\n\n";

exit();
