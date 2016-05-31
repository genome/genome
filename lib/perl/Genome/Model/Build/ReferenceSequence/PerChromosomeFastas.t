#!/usr/bin/env genome-perl

use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
};

use above 'Genome';

use File::Basename qw();
use Genome::Test::Factory::Model::ReferenceSequence;

use Test::More;
use Test::Deep qw(cmp_bag);

my $pkg = 'Genome::Model::Build::ReferenceSequence::PerChromosomeFastas';
use_ok($pkg);

my $test_dir = __FILE__.'.d';

my $reference_fasta_path = File::Spec->join($test_dir,'all_sequences.fa');
my $reference_fasta_absolute_path = File::Spec->rel2abs($reference_fasta_path);

my $reference = Genome::Test::Factory::Model::ReferenceSequence->setup_reference_sequence_build();
my $override = Sub::Override->new(
    'Genome::Model::Build::ReferenceSequence::full_consensus_path',
    sub { return $reference_fasta_absolute_path; }
);

my $result = $pkg->create(
    reference_sequence_build => $reference
);
isa_ok($result, $pkg, 'created per-chromosome-fastas');

my $fasta_list = $result->fasta_list;
ok( scalar(@$fasta_list) == 2, '2 FASTAs created, one for each chromosome');

done_testing();