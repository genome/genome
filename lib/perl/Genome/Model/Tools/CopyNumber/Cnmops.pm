package Genome::Model::Tools::CopyNumber::Cnmops;

use strict;
use warnings;
use Genome;

class Genome::Model::Tools::CopyNumber::Cnmops {
  is => 'Command::V2',
  has_input => [
  tumor_refalign => {
    is => 'Genome::Model::ReferenceAlignment',
    doc => 'Tumor refalign model',
    is_optional => 1,
  },
  normal_refalign => {
    is => 'Genome::Model::ReferenceAlignment',
    doc => 'Normal refalign model',
    is_optional => 1,
  },
  clinseq_model => {
    is => 'Genome::Model::ClinSeq',
    doc => 'Clinseq model(the ref-align models will be extracted from the somatic-exome or somatic-wgs models of this model)',
    is_optional => 1,
  },
  outdir => {
    is => 'FilesystemPath',
    doc => 'Directory to write results',
  },
  roi_bed => {
    is => 'FilesystemPath',
    doc => 'BED file specifying regions to call CNVs on',
    is_optional => 1,
  },
  test => {
    is => 'Boolean',
    doc => 'True for tests, just calls CNVs on some chromosomes',
    default_value => 0,
    is_optional => 1,
  },
  ],
  doc => 'Call CNVs on refalign models(especially Exome) using CnMops',
};

sub help_synopsis {
  return <<EOS
        gmt copy-number cnmops microarray-cnv --outdir=/gscuser/gscuser1/tmp/ --clinseq-model=2887519760
        gmt copy-number cnmops microarray-cnv --outdir=/gscuser/gscuser1/tmp/ --clinseq-model=2887519760 --roi-bed=region_of_interest.bed
        gmt copy-number cnmops microarray-cnv --outdir=/gscuser/gscuser1/tmp/ --tumor-refalign=123456788 --normal-refalign=98765432
EOS
}

sub help_detail {
  return <<EOS
Call CNVs on exome data using CnMops. CnMops has been shown to perform well on WEx data but can be used on WGS data as well. A bed file with the capture regions. A clinseq model with underlying refalign models for normal, tumor  or two refalign models for the normal, tumor have to be provided as input for this tool. A BED file with the regions of interest can be supplied as an optional input, if this is not passed as an input then the ROI file from the refalign build will be used.
EOS
}

sub __errors__ {
  my $self = shift;
  my @errors = $self->SUPER::__errors__(@_);

  unless (-e $self->outdir && -d $self->outdir) {
    push @errors, UR::Object::Tag->create(
      type => 'error',
      properties => ['outdir'],
      desc => "Outdir: " . $self->outdir . " not found or not a directory",
    );
  }
  return @errors;
}

sub get_tumor_normal_refalign_models {
  my $self = shift;
  my $tumor_refalign;
  my $normal_refalign;
  if($self->tumor_refalign and $self->normal_refalign) {
    $tumor_refalign = $self->tumor_refalign;
    $normal_refalign = $self->normal_refalign;
  } elsif($self->clinseq_model) {
    my $base_model;
    if($self->clinseq_model->exome_model) {
      $base_model = $self->clinseq_model->exome_model;
    } else {
      die $self->error_message("Please specify a clinseq model with an exome-sv model");
    }
    if($base_model->tumor_model && $base_model->normal_model) {
      $tumor_refalign = $base_model->tumor_model;
      $normal_refalign = $base_model->normal_model;
    } else {
      die $self->error_message("Please specify a clinseq model with an exome-sv model with refalign builds");
    }
  } else {
    die $self->error_message("Please specify one clinseq model [or] two ref-align models as input!");
  }
  return $tumor_refalign, $normal_refalign;
}

sub get_refalign_bam {
  my $self = shift;
  my $refalign = shift;
  if($refalign->last_succeeded_build->whole_rmdup_bam_file) {
    return $refalign->last_succeeded_build->whole_rmdup_bam_file;
  } else {
    die $self->error_message("Unable to find alignment file for " . $refalign->name);
  }
}

sub get_ROI_bed{
  my $self = shift;
  my $refalign = shift;
  if($self->roi_bed) {
    return $self->roi_bed;
  } elsif($refalign->last_succeeded_build->region_of_interest_set_bed_file) {
    return $refalign->last_succeeded_build->region_of_interest_set_bed_file;
  } else {
    die $self->error_message("Unable to find ROI_set_bed file for " . $refalign->name);
  }
}

sub call_cnmops {
  my $self = shift;
  my $tumor_bam = shift;
  my $normal_bam = shift;
  my $ROI_bed = shift;
  my $outdir = $self->outdir;
  my $cnmops_rscript = __FILE__.".R";
  my $cnmops_rcmd = $cnmops_rscript . " " . $tumor_bam . " " . $normal_bam . " " . $ROI_bed . " " . $outdir ; 
  if($self->test) {
    Genome::Sys->status_message("Running test mode.");
    $cnmops_rcmd = $cnmops_rcmd . " --test";
  } else {
    $cnmops_rcmd = $cnmops_rcmd . " --no-test";
  }
  Genome::Sys->shellcmd(cmd => $cnmops_rcmd);
}

sub execute {
  my $self = shift;
  my $normal_refalign;
  my $tumor_refalign;
  ($tumor_refalign, $normal_refalign) = $self->get_tumor_normal_refalign_models();
  my $tumor_bam = $self->get_refalign_bam($tumor_refalign);
  my $normal_bam = $self->get_refalign_bam($normal_refalign);
  my $ROI_bed = $self->get_ROI_bed($tumor_refalign);
  $self->call_cnmops($tumor_bam, $normal_bam, $ROI_bed);
  $self->status_message("\nCnmops tumor, normal bams are $tumor_bam , $normal_bam");
  $self->status_message("\nROI bed is $ROI_bed");
  return 1;
}

1;
