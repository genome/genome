#!/usr/bin/env genome-perl

use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
};

use above 'Genome';

use Genome::Test::Factory::Build;
use Genome::Test::Factory::Model::CwlPipeline;

use Test::More tests => 11;

my $class = 'Genome::Model::Build::CwlPipeline::InputProcessor';

use_ok($class);

my $model = Genome::Test::Factory::Model::CwlPipeline->setup_object();
my $build = Genome::Test::Factory::Build->setup_object(model_id => $model->id);

for my $n (qw(alpha bravo charlie)) {
    Genome::Model::Build::Input->create(
        name => $n,
        value_id => 'some string value ' . $n,
        value_class_name => 'UR::Value::Text',
        build_id => $build->id,
    );
}

my $i =0;
for my $n (qw(first second third)) {
    Genome::Model::Build::Input->create(
        name => $n,
        value_id => ++$i,
        value_class_name => 'UR::Value::Number',
        build_id => $build->id,
    );
}

for my $n (qw(in out)) {
    Genome::Model::Build::Input->create(
        name => $n,
        value_id => Genome::Sys->create_temp_file_path,
        value_class_name => 'UR::Value::FilePath',
        build_id => $build->id,
    );
}

Genome::Model::Build::Input->create(
    name => 'silly_link_back_to_model',
    value_id => $model->id,
    value_class_name => $model->class,
    build_id => $build->id,
);

{
    package InputProcessor;
    class InputProcessor {
        is => $class,
        has_simple_input => [
            alpha => { input_type => 'Text' },
            bravo => { input_type => 'Text' },
            charlie => { input_type => 'Text' },

            first => { input_type => 'Number' },
            second => { input_type => 'Number' },
            third => { input_type => 'Number' },

            in => { input_type => 'File' },
            out => { input_type => 'File' },
        ],
        has_input => [
            silly_link_back_to_model => { is => 'Genome::Model' },
        ],
    };
}

my $input_processor = InputProcessor->get($build->id);
isa_ok($input_processor, $class);

is($input_processor->build, $build, 'created input processor for build');

my $simple_inputs = $input_processor->simple_inputs;
isa_ok($simple_inputs, 'HASH', 'got back hash of inputs');

my @keys = keys %$simple_inputs;
is(scalar(@keys), 8, 'one key per simple input');

for my $f (qw(in out)) {
    isa_ok($simple_inputs->{$f}, 'HASH');
    is($simple_inputs->{$f}{class}, 'File');
}

ok(!exists $simple_inputs->{silly_link_back_to_model}, 'excluded when not marked simple');
is($input_processor->silly_link_back_to_model, $model, 'object input works');
