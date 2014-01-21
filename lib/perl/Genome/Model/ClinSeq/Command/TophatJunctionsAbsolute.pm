package Genome::Model::ClinSeq::Command::TophatJunctionsAbsolute;
#Written by Malachi Griffith

use strict;
use warnings;
use Genome;
use Genome::Model::ClinSeq::Util qw(:all);

class Genome::Model::ClinSeq::Command::TophatJunctionsAbsolute {
    is => 'Command::V2',
    has_input => [
        build => {
            is => 'Genome::Model::Build::RnaSeq',
            shell_args_position => 1,
            doc => 'RnaSeq build to analyze',
        },
        cancer_annotation_db => {
            is => 'Genome::Db::Tgi::CancerAnnotation',
            example_values => [$Genome::Model::ClinSeq::DEFAULT_CANCER_ANNOTATION_DB_ID],
            doc => 'cancer annotation database'
        },
        outdir => { 
            is => 'FilesystemPath',
            doc => 'Directory where output files will be written', 
        },
    ],
    has_output => [
        junction_topnpercent_file => {
          is => 'FilesystemPath',
          is_optional =>1,
        },
    ],
    doc => 'perform simple processing of tophat junctions output',
};

sub help_synopsis {
    return <<EOS

genome model clin-seq tophat-junctions-absolute --outdir=/tmp/junctions/ 129396808

EOS
}

sub help_detail {
    return <<EOS

This script takes tophat junction results, copies them from rna-seq to clin-seq and performs an post-processing neccessary

EOS
}

sub __errors__ {
  my $self = shift;
  my @errors = $self->SUPER::__errors__;

  #Check output dir
  unless(-e $self->outdir && -d $self->outdir) {
    push @errors,UR::Object::Tag->create(
      type => 'error',
               properties => ['outdir'],
               desc => 'outdir does not exist or is not a directory: ' . $self->outdir,
      );
  }
  return @errors;
}


sub execute {
  my $self = shift;
  my $rnaseq_build = $self->build;
  my $cancer_annotation_db = $self->cancer_annotation_db;
  my $working_dir = $self->outdir;

  $working_dir .= "/" unless ($working_dir =~ /\/$/);
  mkdir ($working_dir) unless (-e $working_dir && -d $working_dir);

  #This step makes a copy of junction results and performs some basic post-processing of the files

  my $rnaseq_build_dir = $rnaseq_build->data_directory;
  my $junctions_dir = $rnaseq_build_dir . "/junctions/";
  if (-e $junctions_dir){
    my $cp_cmd = "cp -r $junctions_dir" . "* $working_dir";
    Genome::Sys->shellcmd(cmd=>$cp_cmd);
  }else{
    $self->error_message("Could not find rna-seq junctions dir for rna-seq build: " . $rnaseq_build->id);
  }

  #Find the file like 'AML103/rnaseq/tumor/tophat_junctions_absolute/NCBI-human.ensembl-67_37l_v2.Junction.GeneExpression.top1percent.tsv' and set as an output to this step
  my $build_id = $rnaseq_build->id;
  my $build_dir = $rnaseq_build->data_directory;
  my $model = $rnaseq_build->model;
  my $ab = $model->annotation_build;
  my $ab_name = $ab->name;
  $ab_name =~ s/\//-/g;
  my $test_path = $working_dir . $ab_name . ".Junction.GeneExpression.top*percent.tsv";
  my $test_result = `ls $test_path`;
  chomp($test_result);
  die $self->error_message("Could not find junctions file in rnaseq build: $build_id ($build_dir)") unless (-e $test_result);

  #Add a 'mapped_gene_name' column to this file and save as an updated version
  my $junction_topnpercent_file = $working_dir . "Junction.GeneExpression.topnpercent.tsv";

  #Get Entrez and Ensembl data for gene name mappings
  my $entrez_ensembl_data = &loadEntrezEnsemblData(-cancer_db => $cancer_annotation_db);

  open(JUNC_IN, "$test_result") || die $self->error_message("Could not open junction file for reading: $test_result");
  open(JUNC_OUT, ">$junction_topnpercent_file") || die $self->error_message("Could not open junction file for writing: $junction_topnpercent_file");
  my $header = 1;
  my %columns;
  while(<JUNC_IN>){
    chomp($_);
    my @line = split("\t", $_);
    if ($header){
      my $p = 0;
      foreach my $column (@line){
        $columns{$column}{position} = $p;
        $p++;
      }
      $header = 0;
      print JUNC_OUT "$_\tmapped_gene_name\n";
      next;
    }
    my $gene_name = $line[$columns{'gene_name'}{position}];
    my $fixed_gene_name = &fixGeneName('-gene'=>$gene_name, '-entrez_ensembl_data'=>$entrez_ensembl_data, '-verbose'=>0);
    print JUNC_OUT "$_\t$fixed_gene_name\n";
  }
  close(JUNC_IN);
  close(JUNC_OUT);

  #Set as an output to this step
  die $self->error_message("Trying to set a file as output but the file does not exist: $junction_topnpercent_file") unless (-e $junction_topnpercent_file);
  $self->junction_topnpercent_file($junction_topnpercent_file);

  return 1;
}

1;
