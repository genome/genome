package Genome::Model::Build::ClinSeq::FileAccessors;

use Genome;
use strict;
use warnings;

class Genome::Model::Build::ClinSeq::FileAccessors {
    is => 'UR::Object',
};

sub case_dir {
  my $self = shift;
  my $case_dir = $self->data_directory . 
    "/" . $self ->common_name;
  return $case_dir; 
}

sub snv_dir {
  my $self = shift;
  my $snv_dir = $self->case_dir . "/snv";
  return $snv_dir;
}

sub cnv_dir {
  my $self = shift;
  my $cnv_dir = $self->case_dir . "/cnv";
  return $cnv_dir;
}

sub microarray_cnv_dir {
  my $self = shift;
  my $microarray_cnv_dir = $self->cnv_dir . "/microarray_cnv";
}

sub exome_snv_dir {
  my $self = shift;
  my $exome_snv_dir = $self->snv_dir . "/exome";
}

sub exome_cnv_dir {
  my $self = shift;
  my $exome_cnv_dir = $self->cnv_dir . "/exome_cnv";
}

sub wgs_snv_dir {
  my $self = shift;
  my $exome_snv_dir = $self->snv_dir . "/wgs";
}

sub wgs_cnv_dir {
  my $self = shift;
  my $wgs_cnv_dir = $self->cnv_dir . "/cnview/CNView_All";
}

sub wgs_exome_snv_dir {
  my $self = shift;
  my $exome_snv_dir = $self->snv_dir . "/wgs_exome";
}

sub wgs_cnvhmm_file {
  my $self = shift;
  my $wgs_cnvhmm_file = $self->wgs_cnv_dir . "/cnaseq.cnvhmm.tsv";
  if(-e $wgs_cnvhmm_file) {
    return $wgs_cnvhmm_file;
  } else {
    die $self->error_message("unable to find $wgs_cnvhmm_file");
  }
}

sub exome_cnvs_file{
  my $self = shift;
  my $exome_cnvs_file = $self->exome_cnv_dir . "/cnmops.cnvs.txt";
  if(-e $exome_cnvs_file) {
    return $exome_cnvs_file;
  } else {
    die $self->error_message("unable to find $exome_cnvs_file");
  }
}

sub microarray_cnvhmm_file {
  my $self = shift;
  my $microarray_cnvhmm_file = $self->microarray_cnv_dir . "/cnvs.diff.cbs.cnvhmm";
  if(-e $microarray_cnvhmm_file) {
    return $microarray_cnvhmm_file
  } else {
    die $self->error_message("unable to find $microarray_cnvhmm_file");
  }
}

1;
