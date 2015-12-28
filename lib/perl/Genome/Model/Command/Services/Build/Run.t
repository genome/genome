#!/usr/bin/env genome-perl

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
}

use strict;
use warnings;

use above 'Genome';

require File::Temp;
use Test::More tests => 8;

use_ok('Genome::Model::Command::Services::Build::Run') or die;

class Genome::Model::Tester { is => 'Genome::ModelDeprecated', };
sub Genome::Model::Tester::_execute_build { return 1; };

my $p = Genome::ProcessingProfile::Tester->create(
    name => 'Tester Test for Testing',
);
ok($p, 'create processing profile');
my $model = Genome::Model->create(
    name => 'Build Run Test',
    processing_profile_id => $p->id,
    subject_type => 'species_name',
    subject_name => 'human'
);
ok($model, 'create model');
my $build = Genome::Model::Build->create(
    id => -7777,
    model => $model
);
ok($build, 'create build');

ok($build->get_or_create_data_directory, 'resolved data dir');

Genome::Report::Email->silent();
ok($build->_initialize_workflow, 'init work flow');
my $run = Genome::Model::Command::Services::Build::Run->execute(
    model_id => $model->id,
    build_id => $build->id,
    inline => 1,
);
ok($run->result, 'run - excute');
is($build->status, 'Succeeded', 'build is succeeded');

done_testing();
