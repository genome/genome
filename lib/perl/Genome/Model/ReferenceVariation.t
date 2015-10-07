#!/usr/bin/env genome-perl

use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
};

use above 'Genome';

use Genome::Test::Factory::ProcessingProfile::ReferenceVariation;
use Genome::Test::Factory::Model::ReferenceSequence;
use Genome::Test::Factory::Sample;

use Test::More tests => 3;

my $pkg = 'Genome::Model::ReferenceVariation';
use_ok($pkg);

my $pp = Genome::Test::Factory::ProcessingProfile::ReferenceVariation->setup_object();
my $reference = Genome::Test::Factory::Model::ReferenceSequence->setup_reference_sequence_build();
my $sample = Genome::Test::Factory::Sample->setup_object;

my $model = $pkg->create(
    processing_profile => $pp,
    name => 'Model/RefereVariation.t test',
    reference_sequence_build => $reference,
    subject => $sample,
);
isa_ok($model, $pkg, 'created model');

is($model->aligner_name, 'speedseq', 'model uses speedseq for alignment');
