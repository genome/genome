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

sub input_dir {
  my $self = shift;
  my $input_dir = $self->case_dir . "/input";
  return $input_dir;
}

sub rnaseq_dir {
  my $self = shift;
  my $rnaseq_dir = $self->case_dir . "/rnaseq";
  return $rnaseq_dir;
}

sub rnaseq_tumor_dir {
  my $self = shift;
  my $rnaseq_tumor_dir = $self->rnaseq_dir . "/tumor";
  return $rnaseq_tumor_dir;
}

sub rnaseq_tumor_cufflinks_dir {
  my $self = shift;
  my $rnaseq_dir = $self->rnaseq_tumor_dir . "/cufflinks_expression_absolute";
  return $rnaseq_dir;
}

sub rnaseq_tumor_cufflinks_genes_dir {
  my $self = shift;
  my $rnaseq_tumor_cufflinks_genes_dir = $self->rnaseq_tumor_cufflinks_dir . "/genes";
  return $rnaseq_tumor_cufflinks_genes_dir;
}

sub rnaseq_tumor_cufflinks_isoforms_dir {
  my $self = shift;
  my $rnaseq_tumor_cufflinks_isoforms_dir = $self->rnaseq_tumor_cufflinks_dir . "/isoforms";
  return $rnaseq_tumor_cufflinks_isoforms_dir;
}

sub rnaseq_tumor_cufflinks_isoforms_merged_dir {
  my $self = shift;
  my $rnaseq_tumor_cufflinks_isoforms_merged_dir = $self->rnaseq_tumor_cufflinks_dir . "/isoforms_merged";
  return $rnaseq_tumor_cufflinks_isoforms_merged_dir;
}

sub rnaseq_tumor_cufflinks_genes_summary_dir {
  my $self = shift;
  my $rnaseq_tumor_cufflinks_genes_summary_dir = $self->rnaseq_tumor_cufflinks_genes_dir . "/summary";
  return $rnaseq_tumor_cufflinks_genes_summary_dir;
}

sub rnaseq_tumor_cufflinks_isoforms_summary_dir {
  my $self = shift;
  my $rnaseq_tumor_cufflinks_isoforms_summary_dir = $self->rnaseq_tumor_cufflinks_isoforms_dir . "/summary";
  return $rnaseq_tumor_cufflinks_isoforms_summary_dir;
}

sub rnaseq_tumor_cufflinks_isoforms_merged_summary_dir {
  my $self = shift;
  my $rnaseq_tumor_cufflinks_isoforms_merged_summary_dir = $self->rnaseq_tumor_cufflinks_isoforms_merged_dir . "/summary";
  return $rnaseq_tumor_cufflinks_isoforms_merged_summary_dir;
}

sub input_summary_dir {
  my $self = shift;
  my $input_summary_dir = $self->input_dir . "/summary";
  return $input_summary_dir;
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

sub exome_snv_summary_dir {
  my $self = shift;
  my $exome_snv_summary_dir =
      $self->exome_snv_dir . "/summary";
  if(-e $exome_snv_summary_dir){
    return $exome_snv_summary_dir;
  } else {
    $self->warning_message("unable to find " .
        $exome_snv_summary_dir);
    return 0;
  }
}

sub wgs_snv_dir {
  my $self = shift;
  my $exome_snv_dir = $self->snv_dir . "/wgs";
}

sub wgs_cnv_dir {
  my $self = shift;
  my $wgs_cnv_dir = $self->cnv_dir . "/wgs_cnv";
}

sub wgs_cnv_cnview_dir {
  my $self = shift;
  my $wgs_cnv_cnview_dir = $self->wgs_cnv_dir . "/cnview/CNView_All";
}

sub wgs_cnv_summary_dir {#CHANGEME
  my $self = shift;
  my $wgs_cnv_summary_dir = $self->cnv_dir . "/wgs_cnv";
}

sub wgs_snv_summary_dir {
  my $self = shift;
  my $wgs_snv_summary_dir =
      $self->wgs_snv_dir . "/summary";
  if(-e $wgs_snv_summary_dir){
    return $wgs_snv_summary_dir;
  } else {
    $self->warning_message("unable to find " .
        $wgs_snv_summary_dir);
    return 0;
  }
}

sub wgs_exome_snv_dir {
  my $self = shift;
  my $exome_snv_dir = $self->snv_dir . "/wgs_exome";
}

sub wgs_exome_snv_summary_dir {
  my $self = shift;
  my $wgs_exome_snv_summary_dir =
      $self->wgs_exome_snv_dir . "/summary";
  if(-e $wgs_exome_snv_summary_dir){
    return $wgs_exome_snv_summary_dir;
  } else {
    $self->warning_message("unable to find " .
        $wgs_exome_snv_summary_dir);
    return 0;
  }
}

sub wgs_cnvhmm_file {
  my $self = shift;
  my $wgs_cnvhmm_file = $self->wgs_cnv_cnview_dir . "/cnaseq.cnvhmm.tsv";
  if(-e $wgs_cnvhmm_file) {
    return $wgs_cnvhmm_file;
  } else {
    $self->warning_message("unable to find $wgs_cnvhmm_file");
    return 0;
  }
}

sub wgs_cnv_wg_plot {
  my $self = shift;
  my $wgs_cnv_wg_plot = $self->wgs_cnv_cnview_dir . "/Both_AllChrs.jpeg";
  if(-e $wgs_cnv_wg_plot) {
    return $wgs_cnv_wg_plot;
  } else {
    $self->warning_message("unable to find $wgs_cnv_wg_plot");
    return 0;
  }
}

sub exome_cnvs_file{
  my $self = shift;
  my $exome_cnvs_file = $self->exome_cnv_dir . "/cnmops.cnvs.txt";
  if(-e $exome_cnvs_file) {
    return $exome_cnvs_file;
  } else {
    $self->warning_message("unable to find $exome_cnvs_file");
    return 0;
  }
}

sub exome_cnv_wg_plot{
  my $self = shift;
  my $exome_cnv_wg_plot = $self->exome_cnv_dir . "/cnmops.segplot_WG.pdf";
  if(-e $exome_cnv_wg_plot) {
    return $exome_cnv_wg_plot;
  } else {
    $self->warning_message("unable to find $exome_cnv_wg_plot");
    return 0;
  }
}

sub microarray_cnvhmm_file {
  my $self = shift;
  my $microarray_cnvhmm_file = $self->microarray_cnv_dir . "/cnvs.diff.cbs.cnvhmm";
  if(-e $microarray_cnvhmm_file) {
    return $microarray_cnvhmm_file
  } else {
    $self->warning_message("unable to find $microarray_cnvhmm_file");
    return 0;
  }
}

sub microarray_cnv_wg_plot {
  my $self = shift;
  my $microarray_cnv_wg_plot = $self->microarray_cnv_dir . "/CNView_All/Gains_AllChrs.jpeg";
  if(-e $microarray_cnv_wg_plot) {
    return $microarray_cnv_wg_plot
  } else {
    $self->warning_message("unable to find $microarray_cnv_wg_plot");
    return 0;
  }
}

sub wgs_exome_snv_tier1_annotated_compact_file {
  my $self = shift;
  my $wgs_exome_snv_tier1_annotated_compact_file =
      $self->wgs_exome_snv_dir . "/snvs.hq.tier1.v1.annotated.compact.tsv";
  if(-e $wgs_exome_snv_tier1_annotated_compact_file ){
    return $wgs_exome_snv_tier1_annotated_compact_file;
  } else {
    $self->warning_message("unable to find " .
        $wgs_exome_snv_tier1_annotated_compact_file);
    return 0;
  }
}

sub wgs_exome_snv_tier1_annotated_compact_catanno_file {
  my $self = shift;
  my $wgs_exome_snv_tier1_annotated_compact_catanno_file =
      $self->wgs_exome_snv_dir . "/snvs.hq.tier1.v1.annotated.compact.catanno.tsv";
  if(-e $wgs_exome_snv_tier1_annotated_compact_catanno_file ){
    return $wgs_exome_snv_tier1_annotated_compact_catanno_file;
  } else {
    $self->warning_message("unable to find " .
        $wgs_exome_snv_tier1_annotated_compact_catanno_file);
    return 0;
  }
}

sub wgs_exome_snv_tier1_annotated_compact_readcounts_file {
  my $self = shift;
  my $wgs_exome_snv_tier1_annotated_compact_readcounts_file =
      $self->wgs_exome_snv_dir . "/snvs.hq.tier1.v1.annotated.compact.readcounts.tsv";
  if(-e $wgs_exome_snv_tier1_annotated_compact_readcounts_file ){
    return $wgs_exome_snv_tier1_annotated_compact_readcounts_file;
  } else {
    $self->warning_message("unable to find " .
        $wgs_exome_snv_tier1_annotated_compact_readcounts_file);
    return 0;
  }
}

sub input_summary_stats_file {
  my $self = shift;
  my $input_summary_stats_file =
      $self->input_summary_dir . "/Stats.tsv";
  if(-e $input_summary_stats_file){
    return $input_summary_stats_file;
  } else {
    $self->warning_message("unable to find " .
        $input_summary_stats_file);
    return 0;
  }
}

sub exome_snv_summary_stats_file {
  my $self = shift;
  my $exome_snv_summary_stats_file =
      $self->exome_snv_summary_dir . "/Stats.tsv";
  if(-e $exome_snv_summary_stats_file){
    return $exome_snv_summary_stats_file;
  } else {
    $self->warning_message("unable to find " .
        $exome_snv_summary_stats_file);
    return 0;
  }
}

sub wgs_exome_snv_summary_stats_file {
  my $self = shift;
  my $wgs_exome_snv_summary_stats_file =
      $self->wgs_exome_snv_summary_dir . "/Stats.tsv";
  if(-e $wgs_exome_snv_summary_stats_file){
    return $wgs_exome_snv_summary_stats_file;
  } else {
    $self->warning_message("unable to find " .
        $wgs_exome_snv_summary_stats_file);
    return 0;
  }
}

sub wgs_snv_summary_stats_file {
  my $self = shift;
  my $wgs_snv_summary_stats_file =
      $self->wgs_snv_summary_dir . "/Stats.tsv";
  if(-e $wgs_snv_summary_stats_file){
    return $wgs_snv_summary_stats_file;
  } else {
    $self->warning_message("unable to find " .
        $wgs_snv_summary_stats_file);
    return 0;
  }
}


sub wgs_cnv_summary_stats_file {
  my $self = shift;
  my $wgs_cnv_summary_stats_file =
      $self->wgs_cnv_summary_dir . "/Stats.tsv";
  if(-e $wgs_cnv_summary_stats_file){
    return $wgs_cnv_summary_stats_file;
  } else {
    $self->warning_message("unable to find " .
        $wgs_cnv_summary_stats_file);
    return 0;
  }
}

sub rnaseq_tumor_cufflinks_genes_stats_file {
  my $self = shift;
  my $rnaseq_tumor_cufflinks_genes_stats_file = $self->rnaseq_tumor_cufflinks_genes_summary_dir .
      "/Stats.tsv";
  if(-e $rnaseq_tumor_cufflinks_genes_stats_file){
    return $rnaseq_tumor_cufflinks_genes_stats_file;
  } else {
    $self->warning_message("unable to find " .
        $rnaseq_tumor_cufflinks_genes_stats_file);
    return 0;
  }
}

sub rnaseq_tumor_cufflinks_isoforms_stats_file {
  my $self = shift;
  my $rnaseq_tumor_cufflinks_isoforms_stats_file = $self->rnaseq_tumor_cufflinks_isoforms_summary_dir .
      "/Stats.tsv";
  if(-e $rnaseq_tumor_cufflinks_isoforms_stats_file){
    return $rnaseq_tumor_cufflinks_isoforms_stats_file;
  } else {
    $self->warning_message("unable to find " .
        $rnaseq_tumor_cufflinks_isoforms_stats_file);
    return 0;
  }
}

sub rnaseq_tumor_cufflinks_isoforms_merged_stats_file {
  my $self = shift;
  my $rnaseq_tumor_cufflinks_isoforms_merged_stats_file = $self->rnaseq_tumor_cufflinks_isoforms_merged_summary_dir .
      "/Stats.tsv";
  if(-e $rnaseq_tumor_cufflinks_isoforms_merged_stats_file){
    return $rnaseq_tumor_cufflinks_isoforms_merged_stats_file;
  } else {
    $self->warning_message("unable to find " .
        $rnaseq_tumor_cufflinks_isoforms_merged_stats_file);
    return 0;
  }
}

1;
