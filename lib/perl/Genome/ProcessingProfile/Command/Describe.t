#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Test::More;

use_ok('Genome::ProcessingProfile::Command::Describe') or die;

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

my $describer = Genome::ProcessingProfile::Command::Describe->create(processing_profiles => [$pp]);
ok($describer, 'Created the describer');
isa_ok($describer, 'Genome::ProcessingProfile::Command::Describe');
ok($describer->execute, 'Executed the describer');

done_testing();
exit;

