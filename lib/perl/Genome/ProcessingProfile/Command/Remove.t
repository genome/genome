#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Test::More;

use_ok('Genome::ProcessingProfile::Command::Remove') or die;

class Genome::Model::Tester {
    is => 'Genome::Model',
};
class Genome::ProcessingProfile::Tester {
    is => 'Genome::ProcessingProfile',
    has_param => [
        param_one =>{ is => 'Text', },
    ],
};
my $pp = Genome::ProcessingProfile::Tester->create(
    name => '__TEST__PP__',
    param_one => 1,
);
ok($pp, "create processing profile") or die;

my $new_name = 'eddie awesome pp for mgc';
my $remover = Genome::ProcessingProfile::Command::Remove->create(
    processing_profile => $pp,
);
ok($remover, 'Created the remover');
$remover->dump_status_messages(1);
ok($remover->execute, 'Executed the remover');

#< BAD >#
# try to remove a pp that has models
$pp = Genome::ProcessingProfile::Tester->create(
    name => '__TEST__PP__',
    param_one => 1,
);
ok($pp, "create processing profile, again") or die;
my $model = Genome::Model::Tester->create(
    name => '__TEST__MODEL__',
    processing_profile => $pp,
    subject => Genome::Sample->create(name => '__TEST__SAMPLE__'),
);
ok($model, "created model") or die;

my $bad = Genome::ProcessingProfile::Command::Remove->create(
    processing_profile => $model->processing_profile,
);
ok($bad, 'Created the remover w/ a processing profile that has a model');
isa_ok($bad, 'Genome::ProcessingProfile::Command::Remove');
ok(!$bad->execute, 'Could not execute the bad remover');

done_testing();
