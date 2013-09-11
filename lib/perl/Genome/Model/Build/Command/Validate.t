#!/usr/bin/env genome-perl
use above "Genome";
use Test::More;

use Genome::Model::TestHelpers qw(
    create_test_sample
    create_test_pp
);

define_test_classes();

my $sample = create_test_sample('test_sample');
my ($build_can, $build_cant ) = create_test_builds("test", $sample);

my $class = "Genome::Model::Build::Command::Validate";
use_ok($class);

my $command_fails = $class->create(builds => [$build_can, $build_cant], builds_can => ["test_property"]);
my $rv;
eval {$rv = $command_fails->execute};
ok (not defined $rv);
ok ($@ =~ m/One or more builds failed validation/);

my $command_succeeds = $class->create(builds => [$build_can, $build_cant], builds_can => ["id"]);
ok ($command_succeeds->execute, "Command succeeded");

done_testing();

sub define_test_classes {
    class Genome::Model::TestCan {
        is => 'Genome::ModelDeprecated',
    };

    class Genome::ProcessingProfile::TestCan {
        is => 'Genome::ProcessingProfile',
    };
    class Genome::Model::Build::TestCan {
        is => 'Genome::Model::Build',
        has => [
            test_property => {
                is => 'String',
                default => "test",
            },
        ],
    };

    class Genome::Model::TestCant {
        is => 'Genome::ModelDeprecated',
    };

    class Genome::ProcessingProfile::TestCant {
        is => 'Genome::ProcessingProfile',
    };
    class Genome::Model::Build::TestCant {
        is => 'Genome::Model::Build',
    };
}

sub create_test_builds {
    my ($base_name, $subject) = @_;

    my @pps, @models, @builds;
    for my $type ("TestCan", "TestCant") {
        my $pp_class = "Genome::ProcessingProfile::$type";
        my $pp = $pp_class->create(
            name => $base_name . $type,
        );
        ok($pp, sprintf('created test pp with name: %s, id: %s',
                $pp->name, $pp->id)) or die;
        push @pps, $pp;

        my $model_class = "Genome::Model::$type";
        my $model = $model_class->create(
            name => $base_name . $type,
            processing_profile => $pp,
            subject => $subject,
        );
        push @models, $model;
        ok($model, sprintf('created test model with name: %s, id: %s',
                $model->name, $model->id)) or die;

        my $build_class = "Genome::Model::Build::$type";
        my $build = $build_class->create(
            model => $model,
        );
        push @builds, $build;
        ok($build, sprintf('created test build id: %s', $build->id)) or die;
    }

    return (@builds);
}
