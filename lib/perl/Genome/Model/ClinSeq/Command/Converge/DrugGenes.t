#!/usr/bin/env perl
use strict;
use warnings;
use above "Genome";
use Test::More tests => 7;
use Cwd;

my $expected_output_dir = $ENV{GENOME_TEST_INPUTS} || '';
ok($expected_output_dir, "found output dir $expected_output_dir");
ok(-e $expected_output_dir, "dir exists");

$expected_output_dir .= '/Genome-Model-ClinSeq-Command-Converge-DrugGenes/2013-02-22/expected-output';
ok($expected_output_dir, "found output dir $expected_output_dir");
ok(-e $expected_output_dir, "dir exists");

my $prefix = UR::Util::used_libs_perl5lib_prefix();
ok($prefix, "using local directory $prefix");

my $cmd_dir = Cwd::abs_path(File::Basename::dirname(__FILE__)) . '/../../original-scripts/';

my $annotation_dir = '/gscmnt/sata132/techd/mgriffit/reference_annotations/';

my $output_dir = Genome::Sys->create_temp_directory();

my $cmd = "PERL5LIB=$prefix:\$PERL5LIB perl $cmd_dir/converge/convergeDrugGenes.pl "
        . " --model_group_id='31779'  --event_types_list='all'  "
        . " --dgidb_subdir_name='drugbank'  --filter_name='antineo'  --outdir=$output_dir" 
        . " --reference_annotations_dir=$annotation_dir  --gene_groups='LUC17'  --verbose=1";

Genome::Sys->shellcmd(cmd => $cmd);
ok(-d $output_dir, "output dir $output_dir");

my @diff = `diff -r --brief $expected_output_dir $output_dir`;
is(scalar(@diff),0, "no differences between expected and actual")
    or do {
        system "mv $output_dir /tmp/last-fail-clinseq-converge-druggenes";
    };

