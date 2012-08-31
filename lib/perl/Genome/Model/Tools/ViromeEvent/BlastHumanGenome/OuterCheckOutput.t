#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More;

use_ok('Genome::Model::Tools::ViromeEvent::BlastHumanGenome::OuterCheckOutput') or die;

my $data_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-ViromeScreening/Titanium17/Titanium17_undecodable';
ok( -d $data_dir, "Test suite dir exists" ) or die;

my $temp_dir = Genome::Sys->create_temp_directory(); #just need some place for log file

my $c = Genome::Model::Tools::ViromeEvent::BlastHumanGenome::OuterCheckOutput->create(
    dir     => $data_dir,
    logfile => $temp_dir.'/log.txt',
    );
ok( $c, "Created blast human genome outer check event" );

ok( $c->execute, "Successfully executed event" );

my $files_for_blast = $c->files_for_blast;
is_deeply( $files_for_blast, [
$ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-ViromeScreening/Titanium17/Titanium17_undecodable/Titanium17_undecodable.fa.cdhit_out.masked.goodSeq_HGblast/Titanium17_undecodable.fa.cdhit_out.masked.goodSeq_file0.fa', ], "Got files for blasting" );

done_testing;

exit;
