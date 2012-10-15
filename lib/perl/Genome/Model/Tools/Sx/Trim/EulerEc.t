#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Test::More;
require File::Compare;

use_ok( 'Genome::Model::Tools::Sx::EulerEc' ) or die;

#input files
my $test_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Sx/';
my $input_file = $test_dir.'/euler.input.fastq';
ok( -s $input_file, "Input file exists" );

#temp dir to run
my $temp_dir = Genome::Sys->create_temp_directory();
ok( -d $temp_dir, "Made temp test dir" );

#output
my $output_file = $temp_dir.'/euler.1_2.output.fastq';

#fail test
my %fail_params1 = (
    input     => [ $input_file.':type=sanger:cnt=2' ],
    output    => [ $output_file.':type=sanger' ],
    min_multi => 10,
);
my $fail_run1 = Genome::Model::Tools::Sx::EulerEc->create( %fail_params1 );
ok( ! $fail_run1->execute, "Failed with kmer_size not set" );
my %fail_params2 = (
    input     => [ $input_file.':type=sanger:cnt=2' ],
    output    => [ $output_file.':type=sanger' ],
    min_multi => 10,
);
my $fail_run2 = Genome::Model::Tools::Sx::EulerEc->create( %fail_params1 );
ok( ! $fail_run2->execute, "Failed with min_multi not set" );

#pass test
my %params = (
    input     => [ $input_file.':type=sanger:cnt=2' ],
    output    => [ $output_file.':type=sanger' ],
    kmer_size => 61,
    min_multi => 10,
);
my $run = Genome::Model::Tools::Sx::EulerEc->create( %params );
ok( $run, 'Created tool' );
ok( $run->execute, 'Executed tool' );

#compare output files
ok( -s $output_file, 'Euler output file created' );
ok( -s $test_dir.'/euler.1_2.output.fastq', 'Example output file exists' );
ok( File::Compare::compare( $output_file,$test_dir.'/euler.1_2.output.fastq' )==0, 'Output file matches example' );

#print $temp_dir."\n";<STDIN>;

done_testing();

exit;
