#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Test::More tests => 12;

BEGIN {
        use_ok('Genome::Model::Tools::Blat::MergeAmplicons');
}

my $tmp_dir = File::Temp::tempdir(CLEANUP => 1);

my $test_data_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Blat-MergeAmplicons';
ok(-d $test_data_dir,'test data dir exists as directory');

my @same_files = glob($test_data_dir .'/same*.txt');
is(scalar(@same_files),3,'Found three same test files');

my @good_merge_files = glob($test_data_dir .'/goodMerge*.txt');
is(scalar(@good_merge_files),3,'Found three good merge test files');

my @bad_merge_files = glob($test_data_dir .'/badMerge*.txt');
is(scalar(@bad_merge_files),2,'Found two bad merge test files');

my $same_file = $tmp_dir .'/same.txt';
my $merge_same_amplicons = Genome::Model::Tools::Blat::MergeAmplicons->create(
                                                                              amplicon_files => \@same_files,
                                                                              output_file => $same_file,
                                                                          );
isa_ok($merge_same_amplicons,'Genome::Model::Tools::Blat::MergeAmplicons');
ok($merge_same_amplicons->execute,'execute command '. $merge_same_amplicons->command_name);

my $good_merge_file = $tmp_dir .'/good_merge.txt';
my $expected_merge_file = $tmp_dir .'/expected_merge_file.txt';
my $expected_cmd = 'cat '. join(' ',@good_merge_files) .' > '. $expected_merge_file;
`$expected_cmd`;

my $good_merge_amplicons = Genome::Model::Tools::Blat::MergeAmplicons->create(
                                                                              amplicon_files => \@good_merge_files,
                                                                              output_file => $good_merge_file,
                                                                          );
isa_ok($merge_same_amplicons,'Genome::Model::Tools::Blat::MergeAmplicons');
ok($good_merge_amplicons->execute,'execute command '. $good_merge_amplicons->command_name);

my $bad_merge_file = $tmp_dir .'/bad_merge.txt';
my $bad_merge_amplicons = Genome::Model::Tools::Blat::MergeAmplicons->create(
                                                                              amplicon_files => \@bad_merge_files,
                                                                              output_file => $bad_merge_file,
                                                                          );
isa_ok($merge_same_amplicons,'Genome::Model::Tools::Blat::MergeAmplicons');
$bad_merge_amplicons->dump_error_messages(0);
$bad_merge_amplicons->queue_error_messages(1);
ok(!$bad_merge_amplicons->execute,'failed execute command '. $bad_merge_amplicons->command_name);
is(scalar(grep { /Different headers found for/ } $bad_merge_amplicons->error_messages),1,'found one error message about different headers');
