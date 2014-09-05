use strict;
use warnings;

use Test::More tests => 1;

use Genome;
use Genome::Test::Factory::DiskAllocation qw();


my $owner = Genome::Sys::User->get(username => Genome::Sys->username);
my $allocation = Genome::Test::Factory::DiskAllocation->setup_object(owner => $owner);

subtest 'finalize adds timeline event' => sub {
    plan tests => 2;
    is(scalar(() = $allocation->timeline_events(name => 'finalized')), 0, 'no finalized event before finalize()');
    $allocation->finalize();
    is(scalar(() = $allocation->timeline_events(name => 'finalized')), 1, 'one finalized event after finalize()');
};
