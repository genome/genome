#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More;

use_ok('Genome::Model::Tools::ViromeEvent::SplitBasedOnBarCode');

my $data_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-ViromeScreening/Titanium17';
ok( -d $data_dir, "Test suite data dir exists" );

my $temp_dir = Genome::Sys->create_temp_directory();

my $sbob = Genome::Model::Tools::ViromeEvent::SplitBasedOnBarCode->create(
    dir         => $data_dir,
    fasta_file  => $data_dir.'/Titanium17_2009_05_05_set0.fna',
    barcode_file=> $data_dir.'/454_Sequencing_log_Titanium_test.txt',
    logfile     => $temp_dir.'/log.txt',
);

ok( $sbob, "Created split-based-on-barcode event" );

done_testing();
