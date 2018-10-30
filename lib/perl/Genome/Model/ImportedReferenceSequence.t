#!/usr/bin/env genome-perl
use strict;
use warnings;
use above "Genome";
use Test::More tests => 6;
use Data::Dumper;

# NOTE this module tests the API accessors, and tests the data behind human build 36
# Creation/building of ImportedReferenceSequence models are tested in Genome/Model/Command/Define/ImportedReferenceSequence.t

BEGIN {
    use_ok("Genome::Model::ImportedReferenceSequence");
}

my $model = Genome::Model::ImportedReferenceSequence->get(name => 'NCBI-human');
isa_ok( $model, 'Genome::Model::ImportedReferenceSequence' );

my $build = $model->build_by_version("36-lite");
my $expected_dir = qr'model_data/2741951221/build101947881$';

like( $build->data_directory(), $expected_dir, 'got the right data directory' );

my $bases_file = $build->get_bases_file(1);
my $expected_bases_file
    = qr'model_data/2741951221/build101947881/1.bases$';
like( $bases_file, $expected_bases_file, 'bases file correct' );

my $test_seq0 = $build->sequence( 1, 1, 10 );
my $expected_seq0 = 'TAACCCTAAC';
is( $test_seq0, $expected_seq0,
    'got expected sequence from start of file' );

my $expected_seq1 = 'NNNNNNNNNNN';
my $test_seq1 = $build->sequence( 1, 247249709, 247249719 );
is( $test_seq1, $expected_seq1,
    'got expected sequence from end of file (via Build...)' );
