#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More;

require File::Compare;

use_ok('Genome::Model::Tools::Mummer::Nucmer') or die;

my $version = 'v1';
my $data_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Mummer/Nucmer-'.$version;
ok ( -d $data_dir, 'Test data dir exists' );

my $test_dir = Genome::Sys->create_temp_directory();
ok ( -d $test_dir, 'Temp test dir created' );

my $query = 'contigs.bases';
my $reference = 'MRSA.MLST.2012-03-09.sequences.fasta';

for my $file ( $query, $reference ) {
    ok ( -s $data_dir."/$file", "Test $file exists" );
    ok ( File::Copy::copy( $data_dir."/$file", $test_dir ), "Copied $file to test dir" );
}

my $tool = Genome::Model::Tools::Mummer::Nucmer->create(
    use_version => '3.22-64',
    prefix      => $test_dir.'/OUT',
    reference   => $test_dir.'/'.$reference,
    query       => $test_dir.'/'.$query,
);
ok( $tool, 'Created nucmer tool' );
ok( $tool->execute, 'Executed nucmer tool' );

my $alignment_out = 'OUT.delta';
ok ( -s $data_dir."/$alignment_out", 'Alignment output file exists' );
ok ( -s $test_dir."/$alignment_out", 'Created new alignment output file' );

ok ( File::Compare::compare( $data_dir."/$alignment_out", $test_dir."/$alignment_out" ) == 1, "Found one line difference in files" );

#print $test_dir."\n";<STDIN>; 

done_testing();
