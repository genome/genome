#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More;

use_ok('Genome::Model::Tools::ViromeEvent::BlastX_NT::InnerCheckOutput') or die;
#check blast dir .. testing on already completed, no-run, blast data so tool won't know if blast db is missing
ok( -s '/gscmnt/sata835/info/medseq/virome/blast_db/nt/nt', "Blast db exists" ) or die;

my $file_to_run = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-ViromeScreening/Titanium17/Titanium17_undecodable/Titanium17_undecodable.BNFiltered_TBLASTX_nt/Titanium17_undecodable.BNFiltered.fa_file0.fa';
ok( -s $file_to_run, "Input file for blast exists" ) or die;

my $done_file = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-ViromeScreening/Titanium17/Titanium17_undecodable/Titanium17_undecodable.BNFiltered_TBLASTX_nt/Titanium17_undecodable.BNFiltered.fa_file0.tblastx.out';
ok( -s $done_file, "Blast output file exists" ) or die; #otherwise will kick off blast which could take a long time

my $temp_dir = Genome::Sys->create_temp_directory();

my $c = Genome::Model::Tools::ViromeEvent::BlastX_NT::InnerCheckOutput->create(
    file_to_run => $file_to_run,
    logfile => $temp_dir.'/log.txt',
);
ok($c, "Created blast-x inner check event") or die;

ok($c->execute, "Successfully executed event");

done_testing;

exit;
