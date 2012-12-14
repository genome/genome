#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::More;
use above "Genome";
require File::Compare;

use_ok( 'Genome::Model::Tools::Newbler::InputsFromFastq' ) or die;

my $version = 'v1';
my $test_suite = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Newbler/InputsFromFastq-'.$version;
ok( -d $test_suite, "Test suite dir exists" ) or die;

my $temp_dir = Genome::Sys->create_temp_directory();
ok( $temp_dir, "Test test dir created" );

my $input_fastq_file = '2869511846-input.fastq';
ok( -s "$test_suite/$input_fastq_file", "Test suite input fastq file exists" ) or die;
symlink( "$test_suite/$input_fastq_file", "$temp_dir/$input_fastq_file" );
ok( -l "$temp_dir/$input_fastq_file", "Input fastq file linked to temp test dir" );

my $create = Genome::Model::Tools::Newbler::InputsFromFastq->create(
    assembly_directory => $temp_dir,
);
ok( $create, "Created tool" );
ok( $create->execute, "Successfully executed tool" );

#compare output files
for my $file ('2869511846-input.fasta', '2869511846-input.fasta.qual' ) {
    ok( -s "$test_suite/consed/edit_dir/$file", "Test suite $file file exists" );
    ok( -s "$temp_dir/consed/edit_dir/$file", "Created $file file" );
    ok( File::Compare::compare("$test_suite/consed/edit_dir/$file","$temp_dir/consed/edit_dir/$file") == 0, "$file files match" );
}

#<STDIN>;

done_testing();
