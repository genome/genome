package Genome::Model::Build::ClinSeq::FileAccessors;

use strict;
use warnings;
use File::Spec;
use Genome;

class Genome::Model::Build::ClinSeq::FileAccessors {
    is => 'UR::Object',
};

sub case_dir {
  my $self = shift;
  my $case_dir = File::Spec->catdir(
    $self->data_directory,
    $self ->common_name);
  return $case_dir;
}

sub snv_dir {
  my $self = shift;
  my $snv_dir = File::Spec->catdir(
    $self->case_dir,
    "snv");
  return $snv_dir;
}

sub indel_dir {
  my $self = shift;
  my $indel_dir = File::Spec->catdir(
    $self->case_dir,
    "indel");
  return $indel_dir;
}

sub snv_indel_report_dir {
  my $self = shift;
  my $bq = shift;
  my $mq = shift;
  my $snv_indel_report_dir = File::Spec->catdir(
    $self->case_dir,
    ("snv_indel_report", "b${bq}_q${mq}"));
  return $snv_indel_report_dir;
}

sub cnv_dir {
  my $self = shift;
  my $cnv_dir = File::Spec->catdir(
    $self->case_dir,
    "cnv");
  return $cnv_dir;
}

sub input_dir {
  my $self = shift;
  my $input_dir = File::Spec->catdir(
    $self->case_dir,
    "input");
  return $input_dir;
}

sub rnaseq_dir {
  my $self = shift;
  my $rnaseq_dir = File::Spec->catdir(
    $self->case_dir,
    "rnaseq");
  return $rnaseq_dir;
}

sub sv_dir {
  my $self = shift;
  my $sv_dir = File::Spec->catdir(
    $self->case_dir,
    "sv");
  return $sv_dir;
}

sub variant_sc_dir {
  my $self = shift;
  my $variant_sc_dir = File::Spec->catdir(
    $self->case_dir,
    "variant_source_callers");
  return $variant_sc_dir;
}

sub mutation_spectrum_dir {
  my $self = shift;
  my $bq = shift;
  my $mq = shift;
  my $mutation_spectrum_dir = File::Spec->catdir(
    $self->case_dir,
    ("mutation-spectrum", "b${bq}_q${mq}"));
  return $mutation_spectrum_dir;
}

sub variant_sc_wgs_dir {
  my $self = shift;
  my $variant_sc_wgs_dir = File::Spec->catdir(
    $self->variant_sc_dir,
    "wgs");
  return $variant_sc_wgs_dir;
}

sub variant_sc_exome_dir {
  my $self = shift;
  my $variant_sc_exome_dir = File::Spec->catdir(
    $self->variant_sc_dir,
    "exome");
  return $variant_sc_exome_dir;
}

sub rnaseq_tumor_dir {
  my $self = shift;
  my $rnaseq_tumor_dir = File::Spec->catdir(
    $self->rnaseq_dir,
    "tumor");
  return $rnaseq_tumor_dir;
}

sub rnaseq_tumor_cufflinks_dir {
  my $self = shift;
  my $rnaseq_dir = File::Spec->catdir(
    $self->rnaseq_tumor_dir,
    "cufflinks_expression_absolute");
  return $rnaseq_dir;
}

sub rnaseq_tumor_cufflinks_genes_dir {
  my $self = shift;
  my $rnaseq_tumor_cufflinks_genes_dir = File::Spec->catdir(
    $self->rnaseq_tumor_cufflinks_dir,
    "genes");
  return $rnaseq_tumor_cufflinks_genes_dir;
}

sub rnaseq_tumor_cufflinks_isoforms_dir {
  my $self = shift;
  my $rnaseq_tumor_cufflinks_isoforms_dir = File::Spec->catdir(
    $self->rnaseq_tumor_cufflinks_dir,
    "isoforms");
  return $rnaseq_tumor_cufflinks_isoforms_dir;
}

sub rnaseq_tumor_cufflinks_isoforms_merged_dir {
  my $self = shift;
  my $rnaseq_tumor_cufflinks_isoforms_merged_dir = File::Spec->catdir(
    $self->rnaseq_tumor_cufflinks_dir,
    "isoforms_merged");
  return $rnaseq_tumor_cufflinks_isoforms_merged_dir;
}

sub rnaseq_tumor_tophat_junctions_absolute_dir {
  my $self = shift;
  my $rnaseq_tumor_tophat_junctions_absolute_dir = File::Spec->catdir(
    $self->rnaseq_tumor_dir,
    "tophat_junctions_absolute");
  return $rnaseq_tumor_tophat_junctions_absolute_dir;
}

sub rnaseq_tumor_cufflinks_genes_summary_dir {
  my $self = shift;
  my $rnaseq_tumor_cufflinks_genes_summary_dir = File::Spec->catdir(
    $self->rnaseq_tumor_cufflinks_genes_dir,
    "summary");
  return $rnaseq_tumor_cufflinks_genes_summary_dir;
}

sub rnaseq_tumor_cufflinks_isoforms_summary_dir {
  my $self = shift;
  my $rnaseq_tumor_cufflinks_isoforms_summary_dir = File::Spec->catdir(
    $self->rnaseq_tumor_cufflinks_isoforms_dir,
    "summary");
  return $rnaseq_tumor_cufflinks_isoforms_summary_dir;
}

sub rnaseq_tumor_cufflinks_isoforms_merged_summary_dir {
  my $self = shift;
  my $rnaseq_tumor_cufflinks_isoforms_merged_summary_dir = File::Spec->catdir(
    $self->rnaseq_tumor_cufflinks_isoforms_merged_dir,
    "summary");
  return $rnaseq_tumor_cufflinks_isoforms_merged_summary_dir;
}

sub rnaseq_tumor_tophat_junctions_absolute_summary_dir {
  my $self = shift;
  my $rnaseq_tumor_tophat_junctions_absolute_summary_dir = File::Spec->catdir(
    $self->rnaseq_tumor_tophat_junctions_absolute_dir,
    "summary");
  return $rnaseq_tumor_tophat_junctions_absolute_summary_dir;
}

sub input_summary_dir {
  my $self = shift;
  my $input_summary_dir = File::Spec->catdir(
    $self->input_dir,
    "summary");
  return $input_summary_dir;
}

sub microarray_cnv_dir {
  my $self = shift;
  my $microarray_cnv_dir = File::Spec->catdir(
    $self->cnv_dir,
    "microarray_cnv");
  return $microarray_cnv_dir;
}

sub exome_snv_dir {
  my $self = shift;
  my $exome_snv_dir = File::Spec->catdir(
    $self->snv_dir,
    "exome");
  return $exome_snv_dir;
}

sub exome_indel_dir {
  my $self = shift;
  my $exome_indel_dir = File::Spec->catdir(
    $self->indel_dir,
    "exome");
  return $exome_indel_dir;
}

sub exome_cnv_dir {
  my $self = shift;
  my $exome_cnv_dir = File::Spec->catdir(
    $self->cnv_dir,
    "exome_cnv");
  return $exome_cnv_dir;
}

sub exome_snv_summary_dir {
  my $self = shift;
  my $exome_snv_summary_dir =
      File::Spec->catdir(
        $self->exome_snv_dir,
        "summary");
  return $exome_snv_summary_dir;
}

sub wgs_snv_dir {
  my $self = shift;
  my $wgs_snv_dir = File::Spec->catdir(
    $self->snv_dir,
    "wgs");
  return $wgs_snv_dir;
}

sub wgs_indel_dir {
  my $self = shift;
  my $wgs_indel_dir = File::Spec->catdir(
    $self->indel_dir,
    "wgs");
  return $wgs_indel_dir;
}

sub wgs_cnv_dir {
  my $self = shift;
  my $wgs_cnv_dir = File::Spec->catdir(
    $self->cnv_dir,
    "wgs_cnv");
  return $wgs_cnv_dir;
}

sub wgs_cnv_cnview_dir {
  my $self = shift;
  my $wgs_cnv_cnview_dir = File::Spec->catdir(
    $self->wgs_cnv_dir,
    qw(cnview CNView_All));
  return $wgs_cnv_cnview_dir;
}

sub wgs_cnv_summary_dir {
  my $self = shift;
  my $wgs_cnv_summary_dir = File::Spec->catdir(
    $self->wgs_cnv_dir,
    "summary");
  return $wgs_cnv_summary_dir;
}

sub wgs_snv_summary_dir {
  my $self = shift;
  my $wgs_snv_summary_dir = File::Spec->catdir(
    $self->wgs_snv_dir,
    "summary");
  return $wgs_snv_summary_dir;
}

sub wgs_exome_indel_dir {
  my $self = shift;
  my $wgs_exome_indel_dir = File::Spec->catdir(
    $self->indel_dir,
    "wgs_exome");
  return $wgs_exome_indel_dir;
}

sub wgs_exome_snv_dir {
  my $self = shift;
  my $wgs_exome_snv_dir = File::Spec->catdir(
    $self->snv_dir,
    "wgs_exome");
  return $wgs_exome_snv_dir;
}

sub wgs_exome_snv_summary_dir {
  my $self = shift;
  my $wgs_exome_snv_summary_dir = File::Spec->catdir(
    $self->wgs_exome_snv_dir,
    "summary");
  return $wgs_exome_snv_summary_dir;
}

sub snv_indel_report_clean_unfiltered_file {
  my $self = shift;
  my $bq = shift;
  my $mq = shift;
  my $snv_indel_report_dir = $self->snv_indel_report_dir($bq, $mq);
  my $snv_indel_report_clean_unfiltered_file = File::Spec->catfile(
    $snv_indel_report_dir,
    $self->common_name . "_final_unfiltered_clean.tsv");
  if(-e $snv_indel_report_clean_unfiltered_file) {
    return $snv_indel_report_clean_unfiltered_file;
  } else {
    $self->warning_message("unable to find " .
      $snv_indel_report_clean_unfiltered_file);
    return;
  }
}

sub snv_indel_report_clean_filtered_file {
  my $self = shift;
  my $bq = shift;
  my $mq = shift;
  my $snv_indel_report_dir = $self->snv_indel_report_dir($bq, $mq);
  my $snv_indel_report_clean_filtered_file = File::Spec->catfile(
    $snv_indel_report_dir,
    $self->common_name . "_final_filtered_clean.tsv");
  if(-e $snv_indel_report_clean_filtered_file) {
    return $snv_indel_report_clean_filtered_file;
  } else {
    $self->warning_message("unable to find " .
      $snv_indel_report_clean_filtered_file);
    return;
  }
}

sub wgs_cnvhmm_file {
  my $self = shift;
  my $wgs_cnvhmm_file = File::Spec->catfile(
    $self->wgs_cnv_cnview_dir,
    "cnaseq.cnvhmm.tsv");
  if(-e $wgs_cnvhmm_file) {
    return $wgs_cnvhmm_file;
  } else {
    $self->warning_message("unable to find $wgs_cnvhmm_file");
    return;
  }
}

sub wgs_cnv_wg_plot {
  my $self = shift;
  my $wgs_cnv_wg_plot = File::Spec->catfile(
    $self->wgs_cnv_cnview_dir,
    "Both_AllChrs.jpeg");
  if(-e $wgs_cnv_wg_plot) {
    return $wgs_cnv_wg_plot;
  } else {
    $self->warning_message("unable to find $wgs_cnv_wg_plot");
    return;
  }
}

sub exome_cnvs_file{
  my $self = shift;
  my $exome_cnvs_file = File::Spec->catfile(
    $self->exome_cnv_dir,
    "cnmops.cnvs.filtered.txt");
  if(-e $exome_cnvs_file) {
    return $exome_cnvs_file;
  } else {
    $self->warning_message("unable to find $exome_cnvs_file");
    return;
  }
}

sub exome_cnv_wg_plot{
  my $self = shift;
  my $exome_cnv_wg_plot = File::Spec->catfile(
    $self->exome_cnv_dir,
    "cnmops.segplot_WG.pdf");
  if(-e $exome_cnv_wg_plot) {
    return $exome_cnv_wg_plot;
  } else {
    $self->warning_message("unable to find $exome_cnv_wg_plot");
    return;
  }
}

sub exome_cnvhmm_file{
  my $self = shift;
  my $exome_cnvhmm_file = File::Spec->catfile(
    $self->exome_cnv_dir,
    "cnmops.cnvhmm");
  if(-e $exome_cnvhmm_file) {
    return $exome_cnvhmm_file;
  } else {
    $self->warning_message("unable to find $exome_cnvhmm_file");
    return;
  }
}

sub microarray_cnvhmm_file {
  my $self = shift;
  my $microarray_cnvhmm_file = File::Spec->catfile(
    $self->microarray_cnv_dir,
    "cnvs.diff.cbs.cnvhmm");
  if(-e $microarray_cnvhmm_file) {
    return $microarray_cnvhmm_file
  } else {
    $self->warning_message("unable to find $microarray_cnvhmm_file");
    return;
  }
}

sub microarray_cnv_wg_plot {
  my $self = shift;
  my $microarray_cnv_wg_plot = File::Spec->catfile(
    $self->microarray_cnv_dir,
    qw(CNView_All Gains_AllChrs.jpeg));
  if(-e $microarray_cnv_wg_plot) {
    return $microarray_cnv_wg_plot
  } else {
    $self->warning_message("unable to find $microarray_cnv_wg_plot");
    return;
  }
}

sub best_cnvhmm_file {
    my $self = shift;
    my $best_cnvhmm_file = $self->wgs_cnvhmm_file;
    unless ($best_cnvhmm_file) {
        $best_cnvhmm_file = $self->exome_cnvhmm_file;
    }
    unless ($best_cnvhmm_file) {
        $best_cnvhmm_file = $self->microarray_cnvhmm_file;
    }
    unless ($best_cnvhmm_file) {
        $self->warning_message("unable to find cnvhmm files");
        return;
    }
    return $best_cnvhmm_file;
}

sub wgs_exome_snv_tier1_annotated_compact_file {
  my $self = shift;
  my $wgs_exome_snv_tier1_annotated_compact_file =
      File::Spec->catfile(
        $self->wgs_exome_snv_dir,
        "snvs.hq.tier1.v1.annotated.compact.tsv");
  if(-e $wgs_exome_snv_tier1_annotated_compact_file ){
    return $wgs_exome_snv_tier1_annotated_compact_file;
  } else {
    $self->warning_message("unable to find " .
        $wgs_exome_snv_tier1_annotated_compact_file);
    return;
  }
}

sub wgs_snv_tier1_annotated_compact_catanno_file {
  my $self = shift;
  my $wgs_snv_tier1_annotated_compact_catanno_file =
      File::Spec->catfile(
        $self->wgs_snv_dir,
        "snvs.hq.tier1.v1.annotated.compact.catanno.tsv");
  if(-e $wgs_snv_tier1_annotated_compact_catanno_file ){
    return $wgs_snv_tier1_annotated_compact_catanno_file;
  } else {
    $self->warning_message("unable to find " .
        $wgs_snv_tier1_annotated_compact_catanno_file);
    return;
  }
}

sub exome_snv_tier1_annotated_compact_catanno_file {
  my $self = shift;
  my $exome_snv_tier1_annotated_compact_catanno_file =
      File::Spec->catfile(
        $self->exome_snv_dir,
        "snvs.hq.tier1.v1.annotated.compact.catanno.tsv");
  if(-e $exome_snv_tier1_annotated_compact_catanno_file ){
    return $exome_snv_tier1_annotated_compact_catanno_file;
  } else {
    $self->warning_message("unable to find " .
        $exome_snv_tier1_annotated_compact_catanno_file);
    return;
  }
}

sub wgs_exome_snv_tier1_annotated_compact_catanno_file {
  my $self = shift;
  my $wgs_exome_snv_tier1_annotated_compact_catanno_file =
      File::Spec->catfile(
        $self->wgs_exome_snv_dir,
        "snvs.hq.tier1.v1.annotated.compact.catanno.tsv");
  if(-e $wgs_exome_snv_tier1_annotated_compact_catanno_file ){
    return $wgs_exome_snv_tier1_annotated_compact_catanno_file;
  } else {
    $self->warning_message("unable to find " .
        $wgs_exome_snv_tier1_annotated_compact_catanno_file);
    return;
  }
}

sub wgs_indel_tier1_annotated_compact_catanno_file {
  my $self = shift;
  my $wgs_indel_tier1_annotated_compact_catanno_file =
      File::Spec->catfile(
        $self->wgs_indel_dir,
        "indels.hq.tier1.v1.annotated.compact.catanno.tsv");
  if(-e $wgs_indel_tier1_annotated_compact_catanno_file ){
    return $wgs_indel_tier1_annotated_compact_catanno_file;
  } else {
    $self->warning_message("unable to find " .
        $wgs_indel_tier1_annotated_compact_catanno_file);
    return;
  }
}

sub exome_indel_tier1_annotated_compact_catanno_file {
  my $self = shift;
  my $exome_indel_tier1_annotated_compact_catanno_file =
      File::Spec->catfile(
        $self->exome_indel_dir,
        "indels.hq.tier1.v1.annotated.compact.catanno.tsv");
  if(-e $exome_indel_tier1_annotated_compact_catanno_file ){
    return $exome_indel_tier1_annotated_compact_catanno_file;
  } else {
    $self->warning_message("unable to find " .
        $exome_indel_tier1_annotated_compact_catanno_file);
    return;
  }
}

sub wgs_exome_indel_tier1_annotated_compact_catanno_file {
  my $self = shift;
  my $wgs_exome_indel_tier1_annotated_compact_catanno_file =
      File::Spec->catfile(
        $self->wgs_exome_indel_dir,
        "indels.hq.tier1.v1.annotated.compact.catanno.tsv");
  if(-e $wgs_exome_indel_tier1_annotated_compact_catanno_file ){
    return $wgs_exome_indel_tier1_annotated_compact_catanno_file;
  } else {
    $self->warning_message("unable to find " .
        $wgs_exome_indel_tier1_annotated_compact_catanno_file);
    return;
  }
}

sub wgs_exome_snv_tier1_annotated_compact_readcounts_file {
  my $self = shift;
  my $wgs_exome_snv_tier1_annotated_compact_readcounts_file =
      File::Spec->catfile(
        $self->wgs_exome_snv_dir,
        "snvs.hq.tier1.v1.annotated.compact.readcounts.tsv");
  if(-e $wgs_exome_snv_tier1_annotated_compact_readcounts_file ){
    return $wgs_exome_snv_tier1_annotated_compact_readcounts_file;
  } else {
    $self->warning_message("unable to find " .
        $wgs_exome_snv_tier1_annotated_compact_readcounts_file);
    return;
  }
}

sub input_summary_stats_file {
  my $self = shift;
  my $input_summary_stats_file =
      File::Spec->catfile(
        $self->input_summary_dir,
        "Stats.tsv");
  if(-e $input_summary_stats_file){
    return $input_summary_stats_file;
  } else {
    $self->warning_message("unable to find " .
        $input_summary_stats_file);
    return;
  }
}

sub exome_snv_summary_stats_file {
  my $self = shift;
  my $exome_snv_summary_stats_file =
      File::Spec->catfile(
        $self->exome_snv_summary_dir,
        "Stats.tsv");
  if(-e $exome_snv_summary_stats_file){
    return $exome_snv_summary_stats_file;
  } else {
    $self->warning_message("unable to find " .
        $exome_snv_summary_stats_file);
    return;
  }
}

sub wgs_exome_snv_summary_stats_file {
  my $self = shift;
  my $wgs_exome_snv_summary_stats_file =
      File::Spec->catfile(
        $self->wgs_exome_snv_summary_dir,
        "Stats.tsv");
  if(-e $wgs_exome_snv_summary_stats_file){
    return $wgs_exome_snv_summary_stats_file;
  } else {
    $self->warning_message("unable to find " .
        $wgs_exome_snv_summary_stats_file);
    return;
  }
}

sub wgs_snv_summary_stats_file {
  my $self = shift;
  my $wgs_snv_summary_stats_file =
      File::Spec->catfile(
        $self->wgs_snv_summary_dir,
        "Stats.tsv");
  if(-e $wgs_snv_summary_stats_file){
    return $wgs_snv_summary_stats_file;
  } else {
    $self->warning_message("unable to find " .
        $wgs_snv_summary_stats_file);
    return;
  }
}


sub wgs_cnv_summary_stats_file {
  my $self = shift;
  my $wgs_cnv_summary_stats_file =
    File::Spec->catfile(
      $self->wgs_cnv_summary_dir,
      "Stats.tsv");
  if(-e $wgs_cnv_summary_stats_file){
    return $wgs_cnv_summary_stats_file;
  } else {
    $self->warning_message("unable to find " .
        $wgs_cnv_summary_stats_file);
    return;
  }
}

sub rnaseq_tumor_cufflinks_genes_stats_file {
  my $self = shift;
  my $rnaseq_tumor_cufflinks_genes_stats_file = File::Spec->catfile(
    $self->rnaseq_tumor_cufflinks_genes_summary_dir,
    "Stats.tsv");
  if(-e $rnaseq_tumor_cufflinks_genes_stats_file){
    return $rnaseq_tumor_cufflinks_genes_stats_file;
  } else {
    $self->warning_message("unable to find " .
        $rnaseq_tumor_cufflinks_genes_stats_file);
    return;
  }
}

sub rnaseq_tumor_cufflinks_isoforms_stats_file {
  my $self = shift;
  my $rnaseq_tumor_cufflinks_isoforms_stats_file = File::Spec->catfile(
    $self->rnaseq_tumor_cufflinks_isoforms_summary_dir,
    "Stats.tsv");
  if(-e $rnaseq_tumor_cufflinks_isoforms_stats_file){
    return $rnaseq_tumor_cufflinks_isoforms_stats_file;
  } else {
    $self->warning_message("unable to find " .
        $rnaseq_tumor_cufflinks_isoforms_stats_file);
    return;
  }
}

sub rnaseq_tumor_cufflinks_isoforms_merged_stats_file {
  my $self = shift;
  my $rnaseq_tumor_cufflinks_isoforms_merged_stats_file = File::Spec->catfile(
    $self->rnaseq_tumor_cufflinks_isoforms_merged_summary_dir,
    "Stats.tsv");
  if(-e $rnaseq_tumor_cufflinks_isoforms_merged_stats_file){
    return $rnaseq_tumor_cufflinks_isoforms_merged_stats_file;
  } else {
    $self->warning_message("unable to find " .
        $rnaseq_tumor_cufflinks_isoforms_merged_stats_file);
    return;
  }
}

sub rnaseq_tumor_tophat_junctions_absolute_stats_file {
  my $self = shift;
  my $rnaseq_tumor_tophat_junctions_absolute_stats_file = File::Spec->catfile(
    $self->rnaseq_tumor_tophat_junctions_absolute_summary_dir,
    "Stats.tsv");
  if(-e $rnaseq_tumor_tophat_junctions_absolute_stats_file){
    return $rnaseq_tumor_tophat_junctions_absolute_stats_file;
  } else {
    $self->warning_message("unable to find " .
        $rnaseq_tumor_tophat_junctions_absolute_stats_file);
    return;
  }
}

sub variant_sc_wgs_stats_file {
  my $self = shift;
  my $variant_sc_wgs_stats_file = File::Spec->catfile(
    $self->variant_sc_wgs_dir,
    "Stats.tsv");
  if(-e $variant_sc_wgs_stats_file){
    return $variant_sc_wgs_stats_file;
  } else {
    $self->warning_message("unable to find " .
        $variant_sc_wgs_stats_file);
    return;
  }
}

sub variant_sc_exome_stats_file {
  my $self = shift;
  my $variant_sc_exome_stats_file = File::Spec->catfile(
    $self->variant_sc_exome_dir,
    "Stats.tsv");
  if(-e $variant_sc_exome_stats_file){
    return $variant_sc_exome_stats_file;
  } else {
    $self->warning_message("unable to find " .
        $variant_sc_exome_stats_file);
    return;
  }
}

sub sv_stats_file {
  my $self = shift;
  my $sv_stats_file = File::Spec->catfile(
    $self->sv_dir,
    "Stats.tsv");
  if(-e $sv_stats_file){
    return $sv_stats_file;
  } else {
    $self->warning_message("unable to find " .
        $sv_stats_file);
    return;
  }
}

sub snv_indel_report_stats_file {
  my $self = shift;
  my $bq = shift;
  my $mq = shift;
  my $sireport_stats_file = File::Spec->catfile(
    $self->snv_indel_report_dir($bq, $mq),
    "Stats.tsv");
  if(-e $sireport_stats_file) {
    return $sireport_stats_file;
  } else {
    $self->warning_message("unable to find " .
        $sireport_stats_file);
    return;
  }
}

sub mutation_spectrum_wgs_summary_file {
  my $self = shift;
  my $bq = shift;
  my $mq = shift;
  my $mutation_spectrum_wgs_summary_file = File::Spec->catfile(
    $self->mutation_spectrum_dir($bq, $mq),
    qw(wgs summarize_mutation_spectrum mutation_spectrum.tsv));
  if(-e $mutation_spectrum_wgs_summary_file ){
    return $mutation_spectrum_wgs_summary_file ;
  } else {
    $self->warning_message("unable to find " .
        $mutation_spectrum_wgs_summary_file);
    return;
  }
}

sub mutation_spectrum_exome_summary_file {
  my $self = shift;
  my $bq = shift;
  my $mq = shift;
  my $mutation_spectrum_exome_summary_file = File::Spec->catfile(
    $self->mutation_spectrum_dir($bq, $mq),
    qw(exome summarize_mutation_spectrum mutation_spectrum.tsv));
  if(-e $mutation_spectrum_exome_summary_file ){
    return $mutation_spectrum_exome_summary_file ;
  } else {
    $self->warning_message("unable to find " .
        $mutation_spectrum_exome_summary_file);
    return;
  }
}

1;
