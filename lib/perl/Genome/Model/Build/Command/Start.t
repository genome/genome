#!/usr/bin/env genome-perl

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;   
};

use strict;
use warnings;

use above 'Genome';
use Test::More;
use Sub::Install qw(reinstall_sub);

Genome::Report::Email->silent();

use_ok('Genome::Model::Build::Command::Start') or die;

class Genome::Model::Tester { is => 'Genome::ModelDeprecated', };
class Genome::Model::Build::Tester { is => 'Genome::Model::Build', };
sub  Genome::Model::Build::Tester::start { return 1; };

my $s = Genome::Sample->create(name => 'TEST-' . __FILE__ . "-$$");
ok($s, "made a test sample");

my $p = Genome::ProcessingProfile::Tester->create(
    name => 'Tester Test for Testing',
);
ok($p, "made a test processing profile");

my $m = Genome::Model::Tester->create(
    processing_profile_id => $p->id,
    subject_class_name => ref($s),
    subject_id => $s->id,
);
ok($m, "made a test model");
my $model_id = $m->id;

reinstall_sub({
    into => 'Genome::Model::Command::Services::Build::Run',
    as => 'execute',
    code => sub {
        return 1;
    },
});

my $cmd = Genome::Model::Build::Command::Start->execute(
    models => [$m],
);
ok($cmd->result, "command believes it succeeded");

done_testing();
