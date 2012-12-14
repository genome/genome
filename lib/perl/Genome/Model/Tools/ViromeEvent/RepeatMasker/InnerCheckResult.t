#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More;

use_ok('Genome::Model::Tools::ViromeEvent::RepeatMasker::InnerCheckResult') or die;

my $file_to_run = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-ViromeScreening/Titanium17/Titanium17_undecodable/Titanium17_undecodable.fa.cdhit_out_RepeatMasker/Titanium17_undecodable.fa.cdhit_out.masked.goodSeq_file0.fa';
ok( -s $file_to_run, "Test data file exists" ) or die;

my $temp_dir = Genome::Sys->create_temp_directory();

my $c = Genome::Model::Tools::ViromeEvent::RepeatMasker::InnerCheckResult->create(
    file_to_run => $file_to_run,
    logfile => $temp_dir.'/log.txt',
    );
ok($c, "Successfully created repeat masker inner check result event");

done_testing();
