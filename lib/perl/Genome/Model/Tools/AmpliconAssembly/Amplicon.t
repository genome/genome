#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Test::More;

use_ok('Genome::Model::Tools::AmpliconAssembly::Amplicon') or die;

my %params = (
    name => 'HMPB-aad13e12',
    directory => $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-AmpliconAssembly/edit_dir',
    reads => [qw/ HMPB-aad13e12.b1 HMPB-aad13e12.b2 HMPB-aad13e12.b3 HMPB-aad13e12.b4 HMPB-aad13e12.g1 HMPB-aad13e12.g2 /],
);
my $amplicon = Genome::Model::Tools::AmpliconAssembly::Amplicon->create(%params);
ok($amplicon, 'create amplicon');

sub invalid_params_for_test_class {
    return (
        directory => 'does_not_exist',
    );
}

for my $attr ( keys %params ) {
    my $method = 'get_'.$attr;
    is_deeply($amplicon->$method, $params{$attr}, "Got $attr");
}

ok($amplicon->get_bioseq, 'Got bioseq');
is($amplicon->get_bioseq_source, 'assembly', 'Got source - assembly');
is($amplicon->was_assembled_successfully, 1, 'Assembled successfully');
is($amplicon->is_bioseq_oriented, 0, 'Not oriented');

my $attempted_reads = $params{reads};

my $assembled_reads = $amplicon->get_assembled_reads;
is_deeply($assembled_reads, $attempted_reads, 'Got source');
is($amplicon->get_assembled_read_count, scalar(@$assembled_reads), 'Got source');
my $read_bioseq = $amplicon->get_bioseq_for_raw_read($attempted_reads->[2]);
is($read_bioseq->id, $attempted_reads->[2], 'Got read bioseq for '.$attempted_reads->[2]);
my $processed_bioseq = $amplicon->get_bioseq_for_processed_read($attempted_reads->[4]);
is($processed_bioseq->id, $attempted_reads->[4], 'Got processed bioseq for '.$attempted_reads->[4]);

my $length = $amplicon->successfully_assembled_length;
is($length, 1150, 'Successfully assembled length');
my $cnt = $amplicon->successfully_assembled_read_count;
is($cnt, 2, 'Successfully assembled read count');
is(
    $amplicon->successfully_assembled_requirements_as_string,
    "length >= $length, reads >= $cnt",
    'Successfully assembled reqs string',
);

done_testing();
exit;

