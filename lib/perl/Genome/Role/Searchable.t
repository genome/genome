#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More tests => 10;

class SearchableThing {
    roles => 'Genome::Role::Searchable',
    has => [ 'prop_a' ],
};

class SearchableThing::View::Solr::Xml {
    is => 'Genome::View::Solr::Xml',
    has_field => [
        type => { is => 'Text', default => 'searchable thing' },
        display_type => { is  => 'Text', default => 'SearchableThing', },
        display_icon_url => { is  => 'Text', default => 'searchable_thing_32', },
        display_url0 => { is => 'Text', default => q(searchable_thing.jpg) },
        default_aspects => { is => 'ARRAY', default => ['prop_a'], },
    ],
};

my $o = SearchableThing->create();
ok($o, 'Created SearchableThing');
my($subject_class, $subject_id) = map { $o->$_ } qw(class id);

ok($o->prop_a('changed'), 'Change attribute on SearchableThing');

ok($o->delete, 'Delete SearchableThing');

my @queues = Genome::Search::Queue->is_loaded();
is(scalar(@queues), 3, '3 Genome::Search::Queue objects created');
foreach my $q ( @queues ) {
    is($q->subject_class, $subject_class, 'subject_class');
    is($q->subject_id, $subject_id, 'subject_id');
}
