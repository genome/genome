#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More tests => 12;

# This test was auto-generated because './Model/ClinSeq/Command/SummarizeTier1SnvSupport.pm'
# had no '.t' file beside it.  Please remove this test if you believe it was
# created unnecessarily.  This is a bare minimum test that just compiles Perl
# and the UR class.
use_ok('Genome::Model::ClinSeq::Command::SummarizeTier1SnvSupport');

#Define the location of expected test results and check that it is valid
my $expected_out = $ENV{GENOME_TEST_INPUTS} . 'Genome-Model-ClinSeq-Command-SummarizeTier1SnvSupport/2013-11-30/';
ok(-d $expected_out, "Directory of expected output exists: $expected_out") or die;

#Get an exome somatic variation build
my $exome_build_id = 138012408;
my $exome_build = Genome::Model::Build->get($exome_build_id);
ok($exome_build, "Got exome somatic variation build from id: $exome_build_id") or die;

#Get a cancer annotation db
my $cancer_annotation_db_name = 'tgi/cancer-annotation/human/build37-20130401.1';
my $cancer_annotation_db = Genome::Db::Tgi::CancerAnnotation->get($cancer_annotation_db_name);
ok($cancer_annotation_db, "Got cancer annotation db from name: $cancer_annotation_db_name") or die;

#Get a clin-seq build
my $clinseq_build_id = "1f51d0be09c94b16878b0bbe8c3709be";
my $clinseq_build = Genome::Model::Build->get($clinseq_build_id);
ok($clinseq_build, "Got a clinseq build from id: $clinseq_build_id") or die;

#Get the clin-seq build dir for this build
my $clinseq_build_dir = $clinseq_build->data_directory;
ok(-d $clinseq_build_dir, "ClinSeq build directory exists") or die;

#Make sure the snv positions file is valid
my $exome_positions_file = $clinseq_build_dir . "/AML109/snv/exome/snvs.hq.tier1.v1.annotated.compact.tsv";
ok (-e $exome_positions_file, "Exome SNV positions file exists") or die;

#Make sure the FPKM file is valid
my $tumor_fpkm_file = $clinseq_build_dir . "/AML109/rnaseq/tumor/cufflinks_expression_absolute/isoforms_merged/isoforms.merged.fpkm.expsort.tsv";
ok (-e $tumor_fpkm_file, "Tumor FPKM file exists") or die;

#Create an output directory
my $output_dir = Genome::Sys->create_temp_directory;
ok (-d $output_dir, "Created temp output dir: $output_dir") or die;

#Run tool as follows:
#genome model clin-seq summarize-tier1-snv-support --cancer-annotation-db='tgi/cancer-annotation/human/build37-20130401.1' --exome-build=138012408  --exome-positions-file=/gscmnt/gc12001/info/model_data/2890221854/build1f51d0be09c94b16878b0bbe8c3709be/AML109/snv/exome/snvs.hq.tier1.v1.annotated.compact.tsv  --tumor-fpkm-file=/gscmnt/gc12001/info/model_data/2890221854/build1f51d0be09c94b16878b0bbe8c3709be/AML109/rnaseq/tumor/cufflinks_expression_absolute/isoforms_merged/isoforms.merged.fpkm.expsort.tsv  --output-dir=/tmp/summarize_tier1_snv_support/
my $cmd = Genome::Model::ClinSeq::Command::SummarizeTier1SnvSupport->create(
  cancer_annotation_db => $cancer_annotation_db,
  exome_build => $exome_build,
  exome_positions_file => $exome_positions_file,
  tumor_fpkm_file => $tumor_fpkm_file,
  output_dir => $output_dir,
);
my $result = $cmd->execute();
is($result, 1, 'Testing for successful execution.  Expecting 1.  Got: ' . $result);

#Since we can not diff the pdf files, at least check for file creation...
my $pdf1 = $output_dir . "/exome/summary/Tumor_VAF_AllDataSources_density.pdf";
ok(-s $pdf1, "Found non-zero PDF file Tumor_VAF_AllDataSources_density.pdf");

#Diff the result, if there are differences, store the new result in /tmp for examination
my $temp_dir = "/tmp/last-summarize-tier1-snv-support/";
my @diff = `diff -x '*.pdf' -x '*.R' -x '.stderr' -x '.stdout' -r $expected_out $output_dir`;
is(scalar(@diff), 0, "only expected differences")
or do {
  for (@diff) { diag($_) }
  Genome::Sys->shellcmd(cmd => "rm -fr $temp_dir");
  Genome::Sys->shellcmd(cmd => "mv $output_dir $temp_dir");
};


