#!/usr/bin/env genome-perl
use strict;
use warnings;

$ENV{UR_DBI_NO_COMMIT} = 1;

use Test::More;
use above "Genome";
use Carp::Always;

{
    package Genome::IsTimelineEvent;

    class Genome::IsTimelineEvent {
        is => 'Genome::Timeline::Event'
    };
}

my $class = 'Genome::Timeline::Event';
use_ok($class);

my $subclass_obj = Genome::IsTimelineEvent->create();
my @properties = keys %{Genome::Timeline::Event->_object_properties};

for my $property (@properties) {
    ok($subclass_obj->can($property), "The preprocessor defines the $property property");
}

done_testing();
