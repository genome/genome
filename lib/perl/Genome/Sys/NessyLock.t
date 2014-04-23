use strict;
use warnings;

use above 'Genome';
use Test::More tests => 2;

use List::Util qw(shuffle);

$ENV{GENOME_NESSY_SERVER} = 'http://nessy.gsc.wustl.edu/';

my $resource_name = 'NessLock.t/' . random_string();
diag 'resource = ' . $resource_name;

my $resource = Genome::Sys::NessyLock->lock(
    resource => $resource_name,
    timeout => 5,
    wait_announce_interval => 10,
);
is($resource, $resource_name, 'locked a Nessy resource');

my $unlocked = Genome::Sys::NessyLock->unlock(resource_lock => $resource);
is($unlocked, 1, 'unlocked a Nessy resource');

sub random_string {
    my @chars = map { (shuffle 'a'..'z', 'A'..'Z', 0..9)[0] } 1..10;
    return join('', @chars);
}
