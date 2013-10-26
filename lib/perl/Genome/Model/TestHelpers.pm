package Genome::Model::TestHelpers;

use strict;
use warnings;

use above 'Genome';
use Test::More;

use Exporter 'import';

our @EXPORT_OK = qw(
    define_test_classes
    create_test_sample
    create_test_pp
    create_test_model
);

sub define_test_classes {
    # Create test subclasses of model and processing profile that can be easily instantiated
    class Genome::Model::Test {
        is => 'Genome::ModelDeprecated',
    };

    class Genome::ProcessingProfile::Test {
        is => 'Genome::ProcessingProfile',
    };
}

sub create_test_sample {
    my ($sample_name) = @_;

    my $sample = Genome::Sample->create(
        name => $sample_name,
    );
    ok($sample, sprintf('created test sample with name: %s, id: %s',
            $sample->name, $sample->id)) or die;
    return $sample;
}

sub create_test_pp {
    my ($pp_name) = @_;

    my $pp = Genome::ProcessingProfile::Test->create(
        name => $pp_name,
    );
    ok($pp, sprintf('created test pp with name: %s, id: %s',
            $pp->name, $pp->id)) or die;
    return $pp;
}

sub create_test_model {
    my ($sample, $pp, $name) = @_;

    my $model = Genome::Model::Test->create(
        subject_id => $sample->id,
        subject_class_name => $sample->class,
        processing_profile_id => $pp->id,
        name => $name
    );
    ok($model, sprintf('created test model with name: %s, id: %s',
            $model->name, $model->id)) or die;
    return $model;
}

1;
