#!/usr/bin/env genome-perl

use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
}

use above 'Genome';
use Test::More;

use Genome::Model::ClinSeq::TestData;

my $class = 'Genome::Model::Build::ClinSeq::FileAccessors';
use_ok($class);

my $data = Genome::Model::ClinSeq::TestData->load();
my $build = Genome::Model::Build->get($data->{CLINSEQ_BUILD});
ok($build, 'generated a test build');

my $data_directory = Genome::Sys->create_temp_directory();
$build->data_directory($data_directory);

my $case_dir = $build->data_directory . "/" . $build->common_name;
Genome::Sys->create_directory($case_dir);
is($build->case_dir, $case_dir, 'found case_dir');

is($build->snv_dir, $case_dir . "/snv", "found snv_dir");
is($build->snv_indel_report_dir(20, 30), $case_dir . "/snv_indel_report/b20_q30",
  "found snv_indel_report dir");
is($build->sv_dir, $case_dir . "/sv", "found sv_dir");
is($build->mutation_spectrum_dir(20, 30), $case_dir . "/mutation-spectrum/b20_q30", "found mutation_spectrum_dir");
is($build->variant_sc_dir, $case_dir . "/variant_source_callers", "found variant_sc_dir");
is($build->variant_sc_wgs_dir, $case_dir . "/variant_source_callers/wgs", "found variant_sc_wgs_dir");
is($build->variant_sc_exome_dir, $case_dir . "/variant_source_callers/exome", "found variant_sc_exome_dir");
is($build->rnaseq_dir, $case_dir . "/rnaseq", "found rnaseq_dir");
is($build->rnaseq_tumor_dir, $case_dir . "/rnaseq/tumor", "found rnaseq_tumor_dir");
is($build->rnaseq_tumor_cufflinks_dir, $case_dir .
      "/rnaseq/tumor/cufflinks_expression_absolute", "found rnaseq_cufflinks_dir");
is($build->rnaseq_tumor_cufflinks_genes_dir, $case_dir .
      "/rnaseq/tumor/cufflinks_expression_absolute/genes", "found rnaseq_cufflinks_genes_dir");
is($build->rnaseq_tumor_tophat_junctions_absolute_dir, $case_dir .
      "/rnaseq/tumor/tophat_junctions_absolute", "found rnaseq_tumor_tophat_junctions_absolute_dir");
is($build->rnaseq_tumor_tophat_junctions_absolute_summary_dir, $case_dir .
      "/rnaseq/tumor/tophat_junctions_absolute/summary", "found rnaseq_tumor_tophat_junctions_absolute_summary_dir");
is($build->rnaseq_tumor_cufflinks_isoforms_dir, $case_dir .
      "/rnaseq/tumor/cufflinks_expression_absolute/isoforms", "found rnaseq_cufflinks_isoforms_dir");
is($build->rnaseq_tumor_cufflinks_isoforms_merged_dir, $case_dir .
      "/rnaseq/tumor/cufflinks_expression_absolute/isoforms_merged", "found rnaseq_cufflinks_isoforms-merged_dir");
is($build->rnaseq_tumor_cufflinks_genes_summary_dir, $case_dir .
      "/rnaseq/tumor/cufflinks_expression_absolute/genes/summary", "found rnaseq_cufflinks_genes_summary_dir");
is($build->rnaseq_tumor_cufflinks_isoforms_summary_dir, $case_dir .
      "/rnaseq/tumor/cufflinks_expression_absolute/isoforms/summary", "found rnaseq_cufflinks_isoforms_summary_dir");
is($build->rnaseq_tumor_cufflinks_isoforms_merged_summary_dir, $case_dir .
      "/rnaseq/tumor/cufflinks_expression_absolute/isoforms_merged/summary", "found rnaseq_cufflinks_isoforms-merged_summary_dir");
is($build->cnv_dir, $case_dir . "/cnv", "found cnv_dir");
is($build->input_dir, $case_dir . "/input", "found input_dir");
is($build->microarray_cnv_dir, $case_dir . "/cnv/microarray_cnv", "found microarray_cnv_dir");
is($build->exome_snv_dir, $case_dir . "/snv/exome", "found exome_snv_dir");
is($build->exome_cnv_dir, $case_dir . "/cnv/exome_cnv", "found exome_cnv_dir");
is($build->wgs_snv_dir, $case_dir . "/snv/wgs", "found wgs_snv_dir");
is($build->wgs_cnv_dir, $case_dir . "/cnv/wgs_cnv", "found wgs_cnv_dir");
is($build->wgs_cnv_summary_dir, $case_dir . "/cnv/wgs_cnv/summary", "found wgs_cnv_summary_dir");
is($build->wgs_cnv_cnview_dir, $case_dir . "/cnv/wgs_cnv/cnview/CNView_All", "found wgs_cnv_cnview_dir");
is($build->wgs_exome_snv_dir, $case_dir . "/snv/wgs_exome", "found wgs_exome_snv_dir");

my @paths_to_create = (
    $case_dir . "/snv_indel_report/b20_q30/" . $build->common_name . "_final_unfiltered_clean.tsv",
    $case_dir . "/snv_indel_report/b20_q30/" . $build->common_name . "_final_filtered_clean.tsv",
    $case_dir . "/cnv/wgs_cnv/cnview/CNView_All/cnaseq.cnvhmm.tsv",
    $case_dir . "/cnv/wgs_cnv/cnview/CNView_All/Both_AllChrs.jpeg",
    $case_dir . "/cnv/exome_cnv/cnmops.cnvs.filtered.txt",
    $case_dir . "/cnv/exome_cnv/cnmops.segplot_WG.pdf",
    $case_dir . "/cnv/exome_cnv/cnmops.cnvhmm",
    $case_dir . "/cnv/microarray_cnv/cnvs.diff.cbs.cnvhmm",
    $case_dir . "/cnv/microarray_cnv/CNView_All/Gains_AllChrs.jpeg",
    $case_dir . "/cnv/wgs_cnv/cnview/CNView_All/cnaseq.cnvhmm.tsv",
    $case_dir . "/snv/wgs_exome/snvs.hq.tier1.v1.annotated.compact.tsv",
    $case_dir . "/snv/wgs/snvs.hq.tier1.v1.annotated.compact.catanno.tsv",
    $case_dir . "/snv/exome/snvs.hq.tier1.v1.annotated.compact.catanno.tsv",
    $case_dir . "/snv/wgs_exome/snvs.hq.tier1.v1.annotated.compact.catanno.tsv",
    $case_dir . "/indel/wgs/indels.hq.tier1.v1.annotated.compact.catanno.tsv",
    $case_dir . "/indel/exome/indels.hq.tier1.v1.annotated.compact.catanno.tsv",
    $case_dir . "/indel/wgs_exome/indels.hq.tier1.v1.annotated.compact.catanno.tsv",
    $case_dir . "/snv/wgs_exome/snvs.hq.tier1.v1.annotated.compact.readcounts.tsv",
    $case_dir . "/snv/exome/summary/Stats.tsv",
    $case_dir . "/snv/exome/snvs.hq.tier1.v1.annotated.compact.tsv.dgidb/expert_antineoplastic.tsv",
    $case_dir . "/indel/exome/indels.hq.tier1.v1.annotated.compact.tsv.dgidb/expert_antineoplastic.tsv",
    $case_dir . "/snv/wgs/snvs.hq.tier1.v1.annotated.compact.tsv.dgidb/expert_antineoplastic.tsv",
    $case_dir . "/indel/wgs/indels.hq.tier1.v1.annotated.compact.tsv.dgidb/expert_antineoplastic.tsv",
    $case_dir . "/cnv/wgs_cnv/cnview/cnv.All_genes.amp.tsv.dgidb/expert_antineoplastic.tsv",
    $case_dir . "/cnv/wgs_cnv/cnview/cnv.All_genes.amp.tsv",
    $case_dir . "/snv/wgs/summary/Stats.tsv",
    $case_dir . "/snv/wgs_exome/summary/Stats.tsv",
    $case_dir . "/input/summary/Stats.tsv",
    $case_dir . "/cnv/wgs_cnv/summary/Stats.tsv",
    $case_dir . "/rnaseq/tumor/cufflinks_expression_absolute/genes/summary/Stats.tsv",
    $case_dir . "/rnaseq/tumor/cufflinks_expression_absolute/isoforms/summary/Stats.tsv",
    $case_dir . "/rnaseq/tumor/cufflinks_expression_absolute/isoforms_merged/summary/Stats.tsv",
    $case_dir . "/rnaseq/tumor/tophat_junctions_absolute/summary/Stats.tsv",
    $case_dir . "/rnaseq/tumor/tophat_junctions_absolute/Junction.GeneExpression.topnpercent.tsv",
    $case_dir . "/rnaseq/tumor/cufflinks_expression_absolute/isoforms_merged/isoforms.merged.fpkm.expsort.top1percent.tsv",
    $case_dir . "/rnaseq/tumor/tophat_junctions_absolute/Junction.GeneExpression.topnpercent.tsv.dgidb/expert_antineoplastic.tsv",
    $case_dir . "/rnaseq/tumor/cufflinks_expression_absolute/isoforms_merged/isoforms.merged.fpkm.expsort.top1percent.tsv.dgidb/expert_antineoplastic.tsv",
    $case_dir . "/sv/Stats.tsv",
    $case_dir . "/variant_source_callers/wgs/Stats.tsv",
    $case_dir . "/variant_source_callers/exome/Stats.tsv",
    $case_dir . "/snv_indel_report/b20_q20/Stats.tsv",
    $case_dir . "/mutation-spectrum/b20_q30/wgs/summarize_mutation_spectrum/mutation_spectrum.tsv",
    $case_dir . "/mutation-spectrum/b20_q30/exome/summarize_mutation_spectrum/mutation_spectrum.tsv",
);

for my $path (@paths_to_create) {
    my (undef, $dir, $file) = File::Spec->splitpath($path);

    Genome::Sys->create_directory($dir);
    Genome::Sys->touch($path);
}


is($build->snv_indel_report_clean_unfiltered_file(20, 30), $case_dir . "/snv_indel_report/b20_q30/" .
   $build->common_name . "_final_unfiltered_clean.tsv",
   "found snv_indel_report_clean_unfiltered_file");
is($build->snv_indel_report_clean_filtered_file(20, 30), $case_dir . "/snv_indel_report/b20_q30/" .
   $build->common_name . "_final_filtered_clean.tsv",
   "found snv_indel_report_clean_filtered_file");
is($build->wgs_cnvhmm_file, $case_dir . "/cnv/wgs_cnv/cnview/CNView_All/cnaseq.cnvhmm.tsv", "found wgs_cnvhmm_file");
is($build->wgs_cnv_wg_plot, $case_dir . "/cnv/wgs_cnv/cnview/CNView_All/Both_AllChrs.jpeg", "found wgs_cnv_wg_plot");
is($build->exome_cnvs_file, $case_dir . "/cnv/exome_cnv/cnmops.cnvs.filtered.txt", "found exome_cnvs_file");
is($build->exome_cnv_wg_plot, $case_dir . "/cnv/exome_cnv/cnmops.segplot_WG.pdf", "found exome_cnvs_wg_plot");
is($build->exome_cnvhmm_file, $case_dir . "/cnv/exome_cnv/cnmops.cnvhmm", "found exome_cnvhmm_file");
is($build->microarray_cnvhmm_file, $case_dir . "/cnv/microarray_cnv/cnvs.diff.cbs.cnvhmm", "found microarray_cnvhmm_file");
is($build->microarray_cnv_wg_plot, $case_dir . "/cnv/microarray_cnv/CNView_All/Gains_AllChrs.jpeg", "found microarray_cnv_wg_plot");
is($build->best_cnvhmm_file, $case_dir . "/cnv/wgs_cnv/cnview/CNView_All/cnaseq.cnvhmm.tsv", "found best_cnvhmm_file");
is($build->wgs_exome_snv_tier1_annotated_compact_file, $case_dir . "/snv/wgs_exome/snvs.hq.tier1.v1.annotated.compact.tsv", "found wgs_exome_snv_tier1_annotated_compact_file");
is($build->wgs_snv_tier1_annotated_compact_catanno_file, $case_dir . "/snv/wgs/snvs.hq.tier1.v1.annotated.compact.catanno.tsv", "found wgs_snv_tier1_annotated_compact_catanno_file");
is($build->exome_snv_tier1_annotated_compact_catanno_file, $case_dir . "/snv/exome/snvs.hq.tier1.v1.annotated.compact.catanno.tsv", "found exome_snv_tier1_annotated_compact_catanno_file");
is($build->wgs_exome_snv_tier1_annotated_compact_catanno_file, $case_dir . "/snv/wgs_exome/snvs.hq.tier1.v1.annotated.compact.catanno.tsv", "found wgs_exome_snv_tier1_annotated_compact_catanno_file");
is($build->wgs_indel_tier1_annotated_compact_catanno_file, $case_dir . "/indel/wgs/indels.hq.tier1.v1.annotated.compact.catanno.tsv", "found wgs_indel_tier1_annotated_compact_catanno_file");
is($build->exome_indel_tier1_annotated_compact_catanno_file, $case_dir . "/indel/exome/indels.hq.tier1.v1.annotated.compact.catanno.tsv", "found exome_indel_tier1_annotated_compact_catanno_file");
is($build->wgs_exome_indel_tier1_annotated_compact_catanno_file, $case_dir . "/indel/wgs_exome/indels.hq.tier1.v1.annotated.compact.catanno.tsv", "found wgs_exome_indel_tier1_annotated_compact_catanno_file");
is($build->wgs_exome_snv_tier1_annotated_compact_readcounts_file, $case_dir . "/snv/wgs_exome/snvs.hq.tier1.v1.annotated.compact.readcounts.tsv", "found wgs_exome_snv_tier1_annotated_compact_readcounts_file");
is($build->exome_snv_summary_dir, $case_dir . "/snv/exome/summary", "found exome_snv_summary_dir");
is($build->wgs_snv_summary_dir, $case_dir . "/snv/wgs/summary", "found wgs_snv_summary_dir");
is($build->wgs_exome_snv_summary_dir, $case_dir . "/snv/wgs_exome/summary", "found wgs_exome_snv_summary_dir");
is($build->exome_snv_summary_stats_file, $case_dir . "/snv/exome/summary/Stats.tsv", "found exome_snv_summary_stats_file");
is($build->exome_snv_dgidb_file, $case_dir . "/snv/exome/snvs.hq.tier1.v1.annotated.compact.tsv.dgidb/expert_antineoplastic.tsv", "found exome_snv_dgidb_file");
is($build->exome_indel_dgidb_file, $case_dir . "/indel/exome/indels.hq.tier1.v1.annotated.compact.tsv.dgidb/expert_antineoplastic.tsv", "found exome_indel_dgidb_file");
is($build->wgs_snv_dgidb_file, $case_dir . "/snv/wgs/snvs.hq.tier1.v1.annotated.compact.tsv.dgidb/expert_antineoplastic.tsv", "found wgs_snv_dgidb_file");
is($build->wgs_indel_dgidb_file, $case_dir . "/indel/wgs/indels.hq.tier1.v1.annotated.compact.tsv.dgidb/expert_antineoplastic.tsv", "found wgs_indel_dgidb_file");
is($build->wgs_cnv_dgidb_file, $case_dir . "/cnv/wgs_cnv/cnview/cnv.All_genes.amp.tsv.dgidb/expert_antineoplastic.tsv", "found wgs_cnv_dgidb_file");
is($build->wgs_cnv_annot_file, $case_dir . "/cnv/wgs_cnv/cnview/cnv.All_genes.amp.tsv", "found wgs_cnv_annot_file");
is($build->wgs_snv_summary_stats_file, $case_dir . "/snv/wgs/summary/Stats.tsv", "found wgs_snv_summary_stats_file");
is($build->wgs_exome_snv_summary_stats_file, $case_dir . "/snv/wgs_exome/summary/Stats.tsv", "found wgs_exome_snv_summary_stats_file");
is($build->input_summary_stats_file, $case_dir . "/input/summary/Stats.tsv", "found input_summary_stats_file");
is($build->wgs_cnv_summary_stats_file, $case_dir . "/cnv/wgs_cnv/summary/Stats.tsv", "found wgs_cnv_summary_stats_file");
is($build->rnaseq_tumor_cufflinks_genes_stats_file, $case_dir .
      "/rnaseq/tumor/cufflinks_expression_absolute/genes/summary/Stats.tsv", "found rnaseq_cufflinks_genes_summary_stats_file");
is($build->rnaseq_tumor_cufflinks_isoforms_stats_file, $case_dir .
      "/rnaseq/tumor/cufflinks_expression_absolute/isoforms/summary/Stats.tsv", "found rnaseq_cufflinks_isoforms_summary_stats_file");
is($build->rnaseq_tumor_cufflinks_isoforms_merged_stats_file, $case_dir .
      "/rnaseq/tumor/cufflinks_expression_absolute/isoforms_merged/summary/Stats.tsv", "found rnaseq_cufflinks_isoforms-merged_summary_stats_file");
is($build->rnaseq_tumor_tophat_junctions_absolute_stats_file, $case_dir .
      "/rnaseq/tumor/tophat_junctions_absolute/summary/Stats.tsv", "found rnaseq_tophat_junctions_absolute_summary_stats_file");
is($build->rnaseq_tumor_tophat_dgidb_dir,
  $case_dir . "/rnaseq/tumor/tophat_junctions_absolute/Junction.GeneExpression.topnpercent.tsv.dgidb", "found rnaseq tophat dgidb dir");
is($build->rnaseq_tumor_cufflinks_dgidb_dir,
  $case_dir . "/rnaseq/tumor/cufflinks_expression_absolute/isoforms_merged/isoforms.merged.fpkm.expsort.top1percent.tsv.dgidb",
  "found cufflinks dgidb dir");
is($build->rnaseq_tumor_tophat_annot_file,
  $case_dir . "/rnaseq/tumor/tophat_junctions_absolute/Junction.GeneExpression.topnpercent.tsv", "found rnaseq tophat annot file");
is($build->rnaseq_tumor_cufflinks_annot_file,
  $case_dir . "/rnaseq/tumor/cufflinks_expression_absolute/isoforms_merged/isoforms.merged.fpkm.expsort.top1percent.tsv",
  "found cufflinks annot file");
is($build->rnaseq_tumor_tophat_dgidb_file,
  $case_dir . "/rnaseq/tumor/tophat_junctions_absolute/Junction.GeneExpression.topnpercent.tsv.dgidb/expert_antineoplastic.tsv", "found tophat dgidb file");
is($build->rnaseq_tumor_cufflinks_dgidb_file,
  $case_dir . "/rnaseq/tumor/cufflinks_expression_absolute/isoforms_merged/isoforms.merged.fpkm.expsort.top1percent.tsv.dgidb/expert_antineoplastic.tsv", "found cufflinks dgidb file");
is($build->sv_stats_file, $case_dir . "/sv/Stats.tsv", "found sv_stats_file");
is($build->variant_sc_wgs_stats_file, $case_dir . "/variant_source_callers/wgs/Stats.tsv", "found variant_sc_wgs_stats_file");
is($build->variant_sc_exome_stats_file, $case_dir . "/variant_source_callers/exome/Stats.tsv", "found variant_sc_exome_stats_file");
is($build->snv_indel_report_stats_file(20, 20), $case_dir . "/snv_indel_report/b20_q20/Stats.tsv", "found snv_indel_report_stats_file");
is($build->mutation_spectrum_wgs_summary_file(20, 30), $case_dir . "/mutation-spectrum/b20_q30/wgs/summarize_mutation_spectrum/mutation_spectrum.tsv",
    "found wgs mut-spec summary file");
is($build->mutation_spectrum_exome_summary_file(20, 30), $case_dir . "/mutation-spectrum/b20_q30/exome/summarize_mutation_spectrum/mutation_spectrum.tsv",
    "found exome mut-spec summary file");

done_testing()
