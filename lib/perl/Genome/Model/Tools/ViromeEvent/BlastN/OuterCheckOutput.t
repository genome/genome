#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More;

use_ok('Genome::Model::Tools::ViromeEvent::BlastN::OuterCheckOutput') or die;

my $data_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-ViromeScreening/Titanium17/Titanium17_undecodable';
ok( -d $data_dir, "Test suite dir exists" ) or die;

my $temp_dir = Genome::Sys->create_temp_directory(); #just need some place for log file

my $c = Genome::Model::Tools::ViromeEvent::BlastN::OuterCheckOutput->create(
    dir     => $data_dir,
    logfile => $temp_dir.'/log.txt',
   );

ok( $c, "Created blast-n outer-check-output event" ) or die;

ok( $c->execute, "Successfully executed event" ) or die;

my $files_for_blast = $c->files_for_blast;
my $expected_file = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-ViromeScreening/Titanium17/Titanium17_undecodable/Titanium17_undecodable.HGfiltered_BLASTN/Titanium17_undecodable.HGfiltered.fa_file0.fa';

is_deeply( $files_for_blast, [ $expected_file, ], "Got files for blasting" );

done_testing();
