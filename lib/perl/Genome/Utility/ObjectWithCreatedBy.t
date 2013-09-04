#!/usr/bin/env genome-perl
use strict;
use warnings;

use Test::More;
use above "Genome";

{
    package Genome::HasCreatedBy;

    class Genome::HasCreatedBy {
        is => 'Genome::Utility::ObjectWithCreatedBy',
        has => [
            dummy_val => { is => 'Text' },
        ]
    };
}


use_ok('Genome::Utility::ObjectWithCreatedBy');

my $obj1 = Genome::HasCreatedBy->create(dummy_val => 6);
is(Genome::Sys->username, $obj1->created_by, 'created_by should be automatically set the first time an object is created');

my $obj2 = Genome::HasCreatedBy->create(created_by => 'turkey', dummy_val => 6);
is('turkey', $obj2->created_by, 'the constructor should not overwrite a passed in value');


done_testing();

