#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More;

use_ok('Genome::Model::Tools::ViromeEvent::BlastHumanGenome::InnerCheckOutput') or die;

my $file_to_run = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-ViromeScreening/Titanium17/Titanium17_undecodable/Titanium17_undecodable.fa.cdhit_out.masked.goodSeq_HGblast/Titanium17_undecodable.fa.cdhit_out.masked.goodSeq_file0.fa';
ok( -s $file_to_run, "Test fasta file exists" ) or die;

my $done_file = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-ViromeScreening/Titanium17/Titanium17_undecodable/Titanium17_undecodable.fa.cdhit_out.masked.goodSeq_HGblast/Titanium17_undecodable.fa.cdhit_out.masked.goodSeq_file0.HGblast.out';
ok( -s $done_file, "Blast completed file exists" ) or die; #otherwise will kick off blast which could take a long time

my $temp_dir = Genome::Sys->create_temp_directory();

my $c = Genome::Model::Tools::ViromeEvent::BlastHumanGenome::InnerCheckOutput->create(
    file_to_run => $file_to_run,
    logfile => $temp_dir.'/foo.txt',
    );
ok( $c, "Created blast human genome event" ) or die;

ok( $c->execute, "Successfully executed event" );

done_testing();
