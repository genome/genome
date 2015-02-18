#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More tests => 5;

use_ok('Genome::Model::Build::ReferenceAlignment');

use Genome::Test::Factory::Build;
use Genome::Test::Factory::InstrumentData::Solexa;
use Genome::Test::Factory::Model::ReferenceAlignment;

my $model = Genome::Test::Factory::Model::ReferenceAlignment->setup_object();
my $other_reference = Genome::Test::Factory::Model::ReferenceAlignment->create_reference_sequence_build();

my $build = Genome::Test::Factory::Build->setup_object(model_id => $model->id);

subtest 'no instrument data does not validate' => sub {
    expect_failure('instrument_data');
};

my $instrument_data = Genome::Test::Factory::InstrumentData::Solexa->setup_object(
    clusters => 5,
    is_paired_end => 1,
);
$model->add_instrument_data($instrument_data);
$build->add_instrument_data($instrument_data);

subtest 'with instrument data does validate' => sub {
    expect_success();
};


my $fl = Genome::FeatureList->__define__(
    name => 'test list for reference alignment validation ' . time(),
    reference => $other_reference,
    format => 'true-BED',
);

$model->region_of_interest_set_name($fl->name);
$build->region_of_interest_set_name($fl->name);

subtest 'roi with different reference does not validate' => sub {
    expect_failure('region_of_interest_set_name');
};

$fl->reference($model->reference_sequence_build);

subtest 'roi with same reference does validate' => sub {
    expect_success();
};


sub expect_failure {
    my $property_name = shift;

    my @tags = $build->validate_for_start();
    is(scalar(@tags), 1, 'found an error in validation') or return;
    my @properties = $tags[0]->properties;
    is($properties[0], $property_name, 'expected property found in tag');
}

sub expect_success {
    my @tags = $build->validate_for_start();
    is(scalar(@tags), 0, 'found no errors in validation')
        or diag(Data::Dumper::Dumper @tags);
}
