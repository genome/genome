#!/usr/bin/env genome-perl
use strict;
use warnings;

use Test::More;
use above "Genome";

{
    package Genome::HasTimestamps;

    class Genome::HasTimestamps {
        is => 'Genome::Utility::ObjectWithTimestamps',
        has => [
            dummy_val => { is => 'Text' },
        ]
    };
}


use_ok('Genome::Utility::ObjectWithTimestamps');

my $inherited_obj = Genome::HasTimestamps->create(dummy_val => 6);
ok($inherited_obj->created_at, 'created_at should be automatically set the first time an object is created');
my $old_val = $inherited_obj->updated_at;
sleep(2);
$inherited_obj->dummy_val(2);
ok($old_val ne $inherited_obj->updated_at, 'updated_at should change when the object changes');

done_testing();

