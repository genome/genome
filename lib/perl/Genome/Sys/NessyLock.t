use strict;
use warnings;

use Genome::Sys::NessyLock;
use Test::More tests => 2;

use List::Util qw(shuffle);

my $resource_name = 'NessLock.t/' . random_string();
diag 'resource = ' . $resource_name;

my $n = Genome::Sys::NessyLock->new(
    url => 'http://nessy.gsc.wustl.edu/',
    is_mandatory => 1,
);

my $resource = $n->lock(
    resource => $resource_name,
    timeout => 5,
    wait_announce_interval => 10,
);
is($resource, $resource_name, 'locked a Nessy resource');

my $unlocked = $n->unlock($resource);
is($unlocked, 1, 'unlocked a Nessy resource');

sub random_string {
    my @chars = map { (shuffle 'a'..'z', 'A'..'Z', 0..9)[0] } 1..10;
    return join('', @chars);
}
