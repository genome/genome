#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More;
use Genome::Test::Factory::Model::SomaticValidation;
use Genome::Test::Factory::Model::ImportedReferenceSequence;
use Genome::Test::Factory::Build;
use Test::MockObject::Extends;

use_ok('Genome::Model::Build::SomaticValidation');

my $sv_model = Genome::Test::Factory::Model::SomaticValidation->setup_object();
my $sv_build = Genome::Test::Factory::Build->setup_object(model_id => $sv_model->id);
my $mock_sv_build = Test::MockObject::Extends->new($sv_build);

my $ref_model = Genome::Test::Factory::Model::ImportedReferenceSequence->setup_object();
my $ref_build = Genome::Test::Factory::Build->setup_object(model_id => $ref_model->id);
my $mock_ref_build = Test::MockObject::Extends->new($ref_build);

my $test_feature = Genome::FeatureList->__define__;
$mock_ref_build->mock('test', sub {return $test_feature;});
$mock_sv_build->mock('reference_sequence_build', sub {return $ref_build});

is($sv_build->get_feature_list_from_reference('test'), $test_feature, "got a feature list from the reference");

done_testing();
