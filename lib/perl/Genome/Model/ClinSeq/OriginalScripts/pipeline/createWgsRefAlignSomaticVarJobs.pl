#!/usr/bin/env genome-perl
#Written by Malachi Griffith
#From a list of patient common names
#Create commands to generate and launch run RefAlign and somatic variation models

use strict;
use warnings;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use Data::Dumper;
use above 'Genome';

#Required
my $common_names = '';

GetOptions ('common_names=s'=>\$common_names);

my $usage=<<INFO;
  Example usage: 
  
  createWgsRefAlignSomaticVarJobs.pl  --common_names='LUC1,LUC2,LUC4,LUC6,LUC7,LUC8,LUC9,LUC10,LUC11,LUC12,LUC13,LUC14,LUC15,LUC16,LUC17,LUC18,LUC20'

INFO

unless ($common_names){
  print GREEN, "\n\n$usage\n\n", RESET;
  exit(1);
}


#Define some model input parameters
my $reference_sequence_build = 106942997; #GRCh37-lite (2869585698)
my $reference_sequence_build_shortname = "Build37/hg19";
my $annotation_reference_build = 106409619; #NCBI-human.combined-annotation (2771411739): NCBI-human.combined-annotation/58_37c_v2 = NCBI-human.ensembl + NCBI-human.genbank
my $annotation_reference_build_shortname =  "58_37c_v2";
my $refalign_processing_profile_name = "Nov 2011 Default Reference Alignment"; #Processing Profile: Nov 2011 Default Reference Alignment (2635769)
my $refalign_processing_profile_id = 2635769;
my $refalign_processing_profile_shortname= "Nov 2011 PP";
my $somaticvariation_processing_profile_name = "Nov 2011 Default Somatic Variation WGS"; #Nov 2011 Default Somatic Variation WGS (2642137)
my $somaticvariation_processing_profile_id = 2642137;
my $somaticvariation_processing_profile_shortname = "Nov 2011 PP";
my $previously_discovered_variations_build = 110108854; #NCBI-human-build37-dbsnp132-whitelist (g1k-human-build37)

my @common_names = split(",", $common_names);

for my $common_name (@common_names) {
  my $new_model_count = 0;
  my %normal_subjects;
  my %tumor_subjects;

  #Get an 'individual object using the patient common name
  print "\n\n\n#$common_name";
  my $individual = Genome::Individual->get(
    common_name => $common_name,
  );
  #Get sample objects associated with the individual object
  my @samples = $individual->samples;
  my $scount = scalar(@samples);
  #print BLUE, "\n\tFound $scount samples", RESET;

  #Get additional info for each sample 
  for my $sample (@samples) {
    #Display basic sample info
    my $sample_name = $sample->name || "UNDEF";
    my $extraction_type = $sample->extraction_type || "UNDEF";
    my $sample_common_name = $sample->common_name || "UNDEF";
    my $tissue_desc = $sample->tissue_desc || "UNDEF";
    my $cell_type = $sample->cell_type || "UNDEF";

    #Skip all but 'genomic dna' samples
    unless ($extraction_type eq "genomic dna"){
      next();
    }
    #Skip samples where the sample common name is undefined
    unless ($sample->common_name){
      next();
    }

    if ($sample_common_name eq "normal"){
      $normal_subjects{$sample_name}{type}="normal";
    }else{
      $tumor_subjects{$sample_name}{type}=$sample_common_name;
    }
    #print MAGENTA, "\n\t\tSAMPLE\tCN: $common_name\tSN: $sample_name\tET: $extraction_type\tSCN: $sample_common_name\tTD: $tissue_desc\tCT: $cell_type", RESET;
  }

  my $normal_count = keys %normal_subjects;
  my $tumor_count = keys %tumor_subjects;

  if ($normal_count > 0 && $tumor_count > 0){

    #Set up the Reference alignment models for tumors then normals

    foreach my $subject_name (sort keys %tumor_subjects){
      $new_model_count++;
      my $sample_type = $tumor_subjects{$subject_name}{type};
      my $model_name = "RefAlign - $common_name - $sample_type - $reference_sequence_build_shortname - $annotation_reference_build_shortname - $refalign_processing_profile_shortname";
      print "\n\n# $new_model_count. Tumor RefAlign Model\t$subject_name\t$common_name\t$refalign_processing_profile_name ($refalign_processing_profile_id)";
      print "\ngenome model define reference-alignment  --model-name='$model_name'  --reference-sequence-build=$reference_sequence_build  --annotation-reference-build=$annotation_reference_build  --subject='$subject_name'  --processing-profile-name='$refalign_processing_profile_name'";
      print "\ngenome model instrument-data assign all-compatible --model-id=''";
      print "\ngenome model build start ''";
    }

    foreach my $subject_name (sort keys %normal_subjects){
      $new_model_count++;
      my $sample_type = $normal_subjects{$subject_name}{type};
      my $model_name = "RefAlign - $common_name - $sample_type - $reference_sequence_build_shortname - $annotation_reference_build_shortname - $refalign_processing_profile_shortname";
      print "\n\n# $new_model_count. Normal RefAlign Model\t$subject_name\t$common_name\t$refalign_processing_profile_name ($refalign_processing_profile_id)";
      print "\ngenome model define reference-alignment  --model-name='$model_name'  --reference-sequence-build=$reference_sequence_build  --annotation-reference-build=$annotation_reference_build  --subject='$subject_name'  --processing-profile-name='$refalign_processing_profile_name'";
      print "\ngenome model instrument-data assign all-compatible --model-id=''";
      print "\ngenome model build start ''";
    }

    #Set up the Somatic variation models for each tumor/normal comparison
    foreach my $normal_subject_name (sort keys %normal_subjects){
      my $normal_sample_type = $normal_subjects{$normal_subject_name}{type};
      foreach my $tumor_subject_name (sort keys %tumor_subjects){
        $new_model_count++;
        my $tumor_sample_type = $tumor_subjects{$tumor_subject_name}{type};
        print "\n\n# $new_model_count. SomaticVariation Model\t$tumor_sample_type vs. $normal_sample_type ($tumor_subject_name vs. $normal_subject_name)\t$common_name\t$somaticvariation_processing_profile_name ($somaticvariation_processing_profile_id)";
        my $model_name = "SomVar - $common_name - $reference_sequence_build_shortname - $annotation_reference_build_shortname - $somaticvariation_processing_profile_shortname";
        print "\ngenome model define somatic-variation  --model-name='$model_name'  --tumor-model=''  --normal-model=''  --processing-profile-name='$somaticvariation_processing_profile_name'  --previously-discovered-variations-build='$previously_discovered_variations_build'  --annotation-build='$annotation_reference_build'";
        print "\ngenome model build start ''";

      }
    }
  }
}

print "\n\n";

exit();
