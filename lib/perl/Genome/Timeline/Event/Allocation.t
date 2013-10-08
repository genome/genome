#!/usr/bin/env genome-perl
use strict;
use warnings;

$ENV{UR_DBI_NO_COMMIT} = 1;

use Test::More;
use above "Genome";
use Carp::Always;

my $class = 'Genome::Timeline::Event::Allocation';

use_ok($class);

for my $event_type (@{Genome::Timeline::Event::AllocationEventType->__meta__->property('id')->valid_values}) {
    ok($class->can($event_type), 'it defines creation methods for each valid event type automatically');
}

ok($class->_properties_to_snapshot(), 'it implements the _properties_to_snapshot sub');

my $test_allocation = Genome::Disk::Allocation->__define__();
my $test_event = Genome::Timeline::Event::Allocation->__define__();

my %prop_hash = $class->_properties_to_snapshot();
while (my ($key, $value) = each %prop_hash) {
    ok($test_allocation->can($key), "$key is a valid property of allocations");
    ok($test_event->can($value), "$value is a valid property of allocation timeline events");
};

done_testing();