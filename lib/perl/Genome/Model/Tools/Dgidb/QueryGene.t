#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use File::Compare;
use Test::More;

use_ok('Genome::Model::Tools::Dgidb::QueryGene');

my $test_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Dgidb-QueryGene/';

my $expected_out = $test_dir.'v1/expected.out';
my $output_file  = Genome::Sys->create_temp_file_path('query_gene.out');

my $cmd =Genome::Model::Tools::Dgidb::QueryGene->create(
    output_file => $output_file,
    genes       => 'FLT,KRAS',
);

ok($cmd, 'command created ok');
ok($cmd->execute, 'command completed successfully');

is(compare($output_file, $expected_out), 0, 'Output file created as expected');

my $expected_outputs = $cmd->output_hash_ref->{matchedTerms};

my $reader = Genome::Utility::IO::SeparatedValueReader->create(
    input     => $output_file,
    separator => "\t",
);

my @outputs;
while (my $data = $reader->next) {
    push @outputs, $data;
}

is_deeply(\@outputs, $expected_outputs, 'Array of hash outputs created as expected.');

done_testing();
