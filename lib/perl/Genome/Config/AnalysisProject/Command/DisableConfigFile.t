#!/usr/bin/env genome-perl

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
}

use strict;
use warnings;

use above 'Genome';

use Test::More;

use_ok('Genome::Config::AnalysisProject::Command::DisableConfigFile') or die;

my $profile_item = Genome::Config::Profile::Item->__define__(
    status => 'active',
);
ok($profile_item, 'define config item');
is($profile_item->status, 'active', "profile item status is 'active'");

my $cmd = Genome::Config::AnalysisProject::Command::DisableConfigFile->execute(profile_item => $profile_item);
ok($cmd->result, 'execute command');
is($profile_item->status, 'disabled', "set profile item status to 'disabled'");

done_testing();
