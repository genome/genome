#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More;

use_ok('Genome::Model::Build::ClinSeq::FileAccessors');

my $m = Genome::Model->get(name => 'apipe-test-clinseq-wer');
ok($m, 'got apipe-test-clinseq-wer');

my $b = $m->last_succeeded_build;
ok($b, 'got a succesful build for apipe-test-clinseq-wer');

my $case_dir = $b->data_directory . "/" . $b->common_name;
is($b->case_dir, $case_dir, 'found case_dir');

is($b->snv_dir, $case_dir . "/snv", "found snv_dir");
is($b->cnv_dir, $case_dir . "/cnv", "found cnv_dir");
is($b->microarray_cnv_dir, $case_dir . "/cnv/microarray_cnv", "found microarray_cnv_dir");
is($b->exome_snv_dir, $case_dir . "/snv/exome", "found exome_snv_dir");
is($b->exome_cnv_dir, $case_dir . "/cnv/exome_cnv", "found exome_cnv_dir");
is($b->wgs_snv_dir, $case_dir . "/snv/wgs", "found wgs_snv_dir");
is($b->wgs_cnv_dir, $case_dir . "/cnv/wgs_cnv/cnview/CNView_All", "found wgs_cnv_dir");
is($b->wgs_exome_snv_dir, $case_dir . "/snv/wgs_exome", "found wgs_exome_snv_dir");
is($b->wgs_cnvhmm_file, $case_dir . "/cnv/wgs_cnv/cnview/CNView_All/cnaseq.cnvhmm.tsv", "found wgs_cnvhmm_file");
is($b->wgs_cnv_wg_plot, $case_dir . "/cnv/wgs_cnv/cnview/CNView_All/Both_AllChrs.jpeg", "found wgs_cnv_wg_plot");
is($b->exome_cnvs_file, $case_dir . "/cnv/exome_cnv/cnmops.cnvs.txt", "found exome_cnvs_file");
is($b->exome_cnv_wg_plot, $case_dir . "/cnv/exome_cnv/cnmops.segplot_WG.pdf", "found exome_cnvs_wg_plot");
is($b->microarray_cnvhmm_file, $case_dir . "/cnv/microarray_cnv/cnvs.diff.cbs.cnvhmm", "found microarray_cnvhmm_file");
is($b->microarray_cnv_wg_plot, $case_dir . "/cnv/microarray_cnv/CNView_All/Gains_AllChrs.jpeg", "found microarray_cnv_wg_plot");
is($b->wgs_exome_snv_tier1_annotated_compact_file, $case_dir . "/snv/wgs_exome/snvs.hq.tier1.v1.annotated.compact.tsv", "found wgs_exome_snv_tier1_annotated_compact_file");
is($b->wgs_exome_snv_tier1_annotated_compact_catanno_file, $case_dir . "/snv/wgs_exome/snvs.hq.tier1.v1.annotated.compact.catanno.tsv", "found wgs_exome_snv_tier1_annotated_compact_catanno_file");
is($b->wgs_exome_snv_tier1_annotated_compact_readcounts_file, $case_dir . "/snv/wgs_exome/snvs.hq.tier1.v1.annotated.compact.readcounts.tsv", "found wgs_exome_snv_tier1_annotated_compact_readcounts_file");

done_testing()
