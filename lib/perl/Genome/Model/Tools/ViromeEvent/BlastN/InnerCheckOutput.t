#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More;

use_ok('Genome::Model::Tools::ViromeEvent::BlastN::InnerCheckOutput') or die;

#check blast dir .. testing on already completed, no-run, blast data so tool won't know if blast db is missing
ok( -s '/gscmnt/sata835/info/medseq/virome/blast_db/nt/nt', "BlastN blast db exists" ) or die;

my $file_to_run = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-ViromeScreening/Titanium17/Titanium17_undecodable/Titanium17_undecodable.HGfiltered_BLASTN/Titanium17_undecodable.HGfiltered.fa_file0.fa';
ok( -s $file_to_run, "Test fasta file exists" ) or die;

my $done_file = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-ViromeScreening/Titanium17/Titanium17_undecodable/Titanium17_undecodable.HGfiltered_BLASTN/Titanium17_undecodable.HGfiltered.fa_file0.blastn.out';
ok( -s $done_file, "Blast out file exists" ) or die; #otherwise will kick off blast which could take a long time

my $temp_dir = Genome::Sys->create_temp_directory();

my $c = Genome::Model::Tools::ViromeEvent::BlastN::InnerCheckOutput->create(
    file_to_run => $file_to_run,
    logfile => $temp_dir.'/log.txt',
    );
ok( $c, "Created nt blastN event" ) or die;

ok( $c->execute, "Successfully executed event" );

done_testing();
