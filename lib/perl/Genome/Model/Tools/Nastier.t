#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Test::More;

if (Genome::Config->arch_os ne 'x86_64') {
    plan skip_all => 'requires 64-bit machine';
}

require File::Compare;

use_ok( 'Genome::Model::Tools::Nastier' );

#test suite dir
my $version = 1;
my $data_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Nastier/v'.$version;
ok( -d $data_dir, 'data dir exists' );

#in/out files
my $input = 'chims.NAST';
my $output = 'chims.NAST.cs';
ok( -s $data_dir.'/'.$input, 'input test file exists' );
ok( -s $data_dir.'/'.$output, 'results file exists' );

#temp test dir
my $temp_dir = Genome::Sys->create_temp_directory();
ok( -d $temp_dir, 'created temp dir' );
ok( File::Copy::copy($data_dir.'/'.$input, $temp_dir), 'copied input fasta file' );

#create/execute tool
my $tool = Genome::Model::Tools::Nastier->create(
    query_FASTA => $temp_dir.'/'.$input,
    output_file => $temp_dir.'/'.$output,
);
ok( $tool, 'created tool' );
ok( $tool->execute, 'executed tool' );

#compare out files
ok( File::Compare::compare($data_dir.'/'.$output, $temp_dir.'/'.$output) == 0, 'Output files match' );

#<STDIN>;

done_testing();
