#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Test::More;

use_ok('Genome::Model::Tools::Sx::PhredEnhancedSeqReader') or die;

is(Genome::Model::Tools::Sx::PhredEnhancedSeqReader->type, 'ephred', 'type is ephred');

my $reader = Genome::Model::Tools::Sx::PhredEnhancedSeqReader->create(file => $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Sx/PhredReaderWriter/in.efasta');
ok($reader, 'create');
my $seq = $reader->read;
ok($seq, 'got seq');
is(
    $seq->{seq}, 
    'GGGGAGGGGAAAAAAAAAAGGGGAAAAAAAAAAAAGGGGaGGGGAAAAAAAAGGGGTTCCTT',
    'sequence matches',
);

done_testing();
