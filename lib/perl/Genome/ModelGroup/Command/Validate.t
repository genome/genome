#!/usr/bin/env genome-perl
use above "Genome";
use Test::More;
use Genome::Model::TestHelpers qw(
    define_test_classes
    create_test_sample
    create_test_pp
    create_test_model
);

define_test_classes();
my $sample = create_test_sample('test_sample');
my $sample2 = create_test_sample('test_sample_2');
my $pp = create_test_pp('test_pp');

my $first = create_test_model($sample, $pp, 'first_test_model');
my $second = create_test_model($sample2, $pp, 'second_test_model');
my $third = create_test_model($sample, $pp, 'third_test_model');

my $class = "Genome::ModelGroup::Command::Validate";
use_ok($class);

my $model_group_dup = Genome::ModelGroup->create(name => "model group with dups", models => [$first,$second,$third]);
my $model_group_nodup = Genome::ModelGroup->create(name => "model group without dups", models => [$first,$second]);

my $command1 = $class->create(model_group => $model_group_dup, unique_subjects => 0);
ok($command1->execute, "Command succeeded");

my $command2 = $class->create(model_group => $model_group_nodup, unique_subjects => 1);
ok($command2->execute, "Command succeeded");

my $command3 = $class->create(model_group => $model_group_dup, unique_subjects => 1);
eval {$command3->execute};
ok($@ =~ m/failed validation/, "Command failed as desired");

done_testing();
