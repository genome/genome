
#!/usr/bin/env genome-perl

use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
}

use above 'Genome';
use Test::More;

my $class = 'Genome::Config::AnalysisProject::Command::Deprecate';
use_ok($class);

my $ap = Genome::Config::AnalysisProject->create(
    name => 'Test Project to Deprecate'
);

is($ap->status, 'Pending', 'initial status should be pending');

my $profile_item = Genome::Config::Profile::Item->create(
    analysis_project => $ap,
    status => 'active',
);  

my $cmd = $class->create(analysis_project => $ap);
ok($cmd->execute(), 'succesful command execution');

is($ap->status, 'Deprecated', 'it should set the status to Deprecated');

is($profile_item->status, 'disabled', 'profile item is disabled');     

my $fail_cmd = $class->create(analysis_project => $ap);
ok(!$fail_cmd->execute(), 'fail command execution');

done_testing();

1;
