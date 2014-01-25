#!/usr/bin/env genome-perl
#Written by Malachi Griffith

#This script gather basic stats for a series of optional builds and write out summary results to a working directory
#Optional builds:
#WGS somatic variation build
#Exome somatic variation build
#RNA-seq build (for tumor)
#RNA-seq build (for normal)

#NOTE:  This is distinct from the kinds of things that will be done in Summarize.pm!
#That script works on a completed ClinSeq model/build and tell you all about that
#This script will be run as *part* of the ClinSeq run and will produce a summary of the models/builds/processing profiles etc. that were used


#Load modules
use strict;
use warnings;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use Data::Dumper;
use above 'Genome'; 
use Genome::Model::ClinSeq::Util qw(:all);

#Input parameters
my $wgs_som_var_build_id = '';
my $exome_som_var_build_id = '';
my $tumor_rna_seq_build_id = '';
my $normal_rna_seq_build_id = '';
my $working_dir = '';
my $verbose = 0;
my $clean = 0;

#Example data:
#WGS somatic var   model id = 2880644349  |  build id = 115623342
#Exome somatic var model id = 2880732183  |  build id = 115672230
#Tumor rna-seq     model id = 2880794563  |  build id = 115909755

GetOptions ('tumor_rna_seq_build_id=s'=>\$tumor_rna_seq_build_id, 'normal_rna_seq_build_id=s'=>\$normal_rna_seq_build_id,
	    'wgs_som_var_build_id=s'=>\$wgs_som_var_build_id, 'exome_som_var_build_id=s'=>\$exome_som_var_build_id, 
 	    'working_dir=s'=>\$working_dir, 'verbose=i'=>\$verbose, 'clean=i'=>\$clean);

my $usage=<<INFO;

  Example usage: 
  
  clinseq.pl  --wgs_som_var_build_id='115623342'  --exome_som_var_build_id='115672230'  --tumor_rna_seq_build_id='115909755'  --working_dir=/gscmnt/sata132/techd/mgriffit/hgs/
  
  Intro:
  This script attempts to automate the process of running the 'clinseq' pipeline

  Details:
  --wgs_som_var_build_id          Whole genome sequence (WGS) somatic variation build ID
  --exome_som_var_build_id        Exome capture sequence somatic variation build ID
  --tumor_rna_seq_build_id        RNA-seq build id for the tumor sample
  --normal_rna_seq_build_id       RNA-seq build id for the normal sample
  --working_dir                   Directory where a patient subdir will be created
  --common_name                   Patient's common name (will be used for the name of a results dir and labeling purposes)
  --verbose                       To display more output, set to 1
  --clean                         To clobber the top dir and create everything from scratch, set to 1

INFO

unless (($wgs_som_var_build_id || $exome_som_var_build_id || $tumor_rna_seq_build_id || $normal_rna_seq_build_id) && $working_dir){
  print GREEN, "$usage", RESET;
  exit();
}

#Set flags for each datatype
my ($wgs, $exome, $tumor_rnaseq, $normal_rnaseq) = (0,0,0,0);
if ($wgs_som_var_build_id){$wgs=1;}
if ($exome_som_var_build_id){$exome=1;}
if ($tumor_rna_seq_build_id){$tumor_rnaseq=1;}
if ($normal_rna_seq_build_id){$normal_rnaseq=1;}

#Check the working dir
$working_dir = &checkDir('-dir'=>$working_dir, '-clear'=>"no");

#Get the builds
&getBuilds();

my $builds = &getBuilds('-wgs_som_var_build_id'=>$wgs_som_var_build_id, '-exome_som_var_build_id'=>$exome_som_var_build_id, '-tumor_rna_seq_build_id'=>$tumor_rna_seq_build_id, '-normal_rna_seq_build_id'=>$normal_rna_seq_build_id);

print Dumper $builds;


exit();


########################################################################################################################################################
#Get basic build info from the specified build IDs (e.g. BAM paths, build dirs, reference annotations, etc)                                            #
########################################################################################################################################################
sub getBuilds{
  my %args = @_;
  my $wgs_build_id = $args{'-wgs_som_var_build_id'};
  my $exome_build_id = $args{'-exome_som_var_build_id'};
  my $tumor_rna_build_id = $args{'-tumor_rna_seq_build_id'};
  my $normal_rna_build_id = $args{'-normal_rna_seq_build_id'};
  my %b;

  my ($wgs_datadir, $exome_datadir, $tumor_rna_datadir, $normal_rna_datadir) = ('', '', '', '');
  if ($wgs_build_id){
    my $wgs_build = Genome::Model::Build->get($wgs_build_id);
    if ($wgs_build){
        #... /genome/lib/perl/Genome/Model/Build/SomaticVariation.pm
        $b{wgs}{build_dir} =  $wgs_build->data_directory ."/";
        $b{wgs}{normal_bam} = $wgs_build->normal_bam;
        $b{wgs}{tumor_bam} = $wgs_build->tumor_bam;
        $b{wgs}{model_id} = $wgs_build->model->id;
        $b{wgs}{common_name} = $wgs_build->subject->patient->common_name;

    }else{
      print RED, "\n\nA WGS build ID was specified, but a build object could not be found!\n\n", RESET;
      exit();
    }
  }

  if ($exome_build_id){
    my $exome_build = Genome::Model::Build->get($exome_build_id);
    if ($exome_build){
        #... /genome/lib/perl/Genome/Model/Build/SomaticVariation.pm
        $b{exome}{build_dir} =  $exome_build->data_directory ."/";
        $b{exome}{normal_bam} = $exome_build->normal_bam;
        $b{exome}{tumor_bam} = $exome_build->tumor_bam;
        $b{exome}{model_id} = $exome_build->model->id;
        $b{exome}{common_name} = $exome_build->subject->patient->common_name;

    }else{
      print RED, "\n\nA WGS build ID was specified, but a build object could not be found!\n\n", RESET;
      exit();
    }
  }

  if ($tumor_rna_build_id){
    my $rna_build = Genome::Model::Build->get($tumor_rna_build_id);
    if ($rna_build){
      $b{tumor_rnaseq}{build_dir} = $rna_build->data_directory ."/";
      my $alignment_result = $rna_build->alignment_result;
      $b{tumor_rnaseq}{bam} = $alignment_result->bam_file;
      $b{tumor_rnaseq}{model_id} = $rna_build->model->id;
      $b{tumor_rnaseq}{common_name} = $rna_build->subject->patient->common_name;

    }else{
      print RED, "\n\nA tumor RNA-seq build ID was specified, but a build object could not be found!\n\n", RESET;
      exit();
    }
  }

  if ($normal_rna_build_id){
    my $rna_build = Genome::Model::Build->get($normal_rna_build_id);
    if ($rna_build){
      $b{normal_rnaseq}{build_dir} = $rna_build->data_directory ."/";
      my $alignment_result = $rna_build->alignment_result;
      $b{normal_rnaseq}{bam} = $alignment_result->bam_file;
    }else{
      print RED, "\n\nA normal RNA-seq build ID was specified, but a build object could not be found!\n\n", RESET;
      exit();
    }
  }
 



  return(\%b);
}





