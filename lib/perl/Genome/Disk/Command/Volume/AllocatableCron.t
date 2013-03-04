use strict;
use warnings;

use Test::More;
use above 'Genome';

BEGIN {
    use_ok('Genome::Disk::Volume');
    use_ok('Genome::Disk::Command::Volume::AllocatableCron');
};

my @volume_methods = qw(
    can_allocate
    is_over_soft_limit
    is_used_over_soft_limit
    is_allocated_over_soft_limit
);
for my $m (@volume_methods) {
    ok(Genome::Disk::Volume->can($m), "Genome::Disk::Volume has '$m' method");
}

done_testing();
