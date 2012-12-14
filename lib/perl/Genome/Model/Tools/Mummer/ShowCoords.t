#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More;

require File::Compare;

use_ok('Genome::Model::Tools::Mummer::ShowCoords') or die;

my $version = 'v1';
my $data_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Mummer/ShowCoords-'.$version;
ok( -d $data_dir, "Test data dir exists" );

my $test_dir = Genome::Sys->create_temp_directory();
ok( -d $test_dir, "Temp test dir created" );

my $input_file = 'OUT.delta';
my $output_file = 'alignments.txt';
for my $file ( $input_file, $output_file ) {
    ok( -s $data_dir."/$file", "Test $file exists" );
}

my $tool = Genome::Model::Tools::Mummer::ShowCoords->create(
    input_delta_file => $data_dir."/$input_file",
    output_file => $test_dir."/$output_file",
    r => 1,
    c => 1,
    l => 1,
    T => 1,
);
ok( $tool, 'Created tool' );
ok( $tool->execute, 'Successfully executed tool' );

ok( File::Compare::compare( $data_dir."/$output_file",$test_dir."/$output_file" ) == 0, 'Output files match' );

#<STDIN>;
done_testing();
