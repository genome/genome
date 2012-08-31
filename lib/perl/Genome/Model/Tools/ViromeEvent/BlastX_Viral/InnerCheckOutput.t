#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More;

use_ok('Genome::Model::Tools::ViromeEvent::BlastX_Viral::InnerCheckOutput') or die;
#check blast dir .. testing on already completed, no-run, blast data so tool won't know if blast db is missing
ok( -s '/gscmnt/sata835/info/medseq/virome/blast_db/viral/viral.genomic.fna', "Blast db exists" );

my $file_for_blast = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-ViromeScreening/Titanium17/Titanium17_undecodable/Titanium17_undecodable.TBXNTFiltered_TBLASTX_ViralGenome/Titanium17_undecodable.TBXNTFiltered.fa_file0.fa';
ok( -s $file_for_blast, "File for blast exists" ) or die;
my $done_file = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-ViromeScreening/Titanium17/Titanium17_undecodable/Titanium17_undecodable.TBXNTFiltered_TBLASTX_ViralGenome/Titanium17_undecodable.TBXNTFiltered.fa_file0.tblastx_ViralGenome.out';
ok( -s $done_file, "Blast output file exists" ) or die; #otherwise will kick off blast which could take a long time

my $temp_dir = Genome::Sys->create_temp_directory();

my $c = Genome::Model::Tools::ViromeEvent::BlastX_Viral::InnerCheckOutput->create(
    file_to_run => $file_for_blast,
    logfile => $temp_dir.'/log.txt',
    );
ok( $c, "Created blastx viral inner check output event" ) or die;

ok( $c->execute, "Successfully executed event" );

done_testing();;

exit;
