#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Test::More;

use_ok('Genome::Report::Dataset');

my $dataset = Genome::Report::Dataset->create(
    name => 'stats',
    row_name => 'stat',
    headers => [qw/ attempted assembled assembled-percent /],
    rows => [ [qw/ 4 5 80.0% /] ],
    attributes => { good => 'mizzou', bad => 0 },
);
ok($dataset, 'create dataset');

ok($dataset->to_xml_element, 'XML element');
ok($dataset->to_xml_string, 'XML string');

is(
    $dataset->to_separated_value_string(separator => '|'), 
    "attempted|assembled|assembled-percent\n4|5|80.0\%\n",
    'SVS string',
);
is(
    $dataset->to_separated_value_string(separator => '|', include_headers => 0), 
    "4|5|80.0\%\n",
    'SVS string',
);

is(
    $dataset->get_attribute('good'), 
    'mizzou',
    'get_attribute good => mizzou',
);
ok(
    $dataset->set_attribute('bad', 'kansas'), 
    'set_attribute bad => kansas',
);

my $from_element_ds = Genome::Report::Dataset->create_from_xml_element($dataset->to_xml_element);
ok($from_element_ds, 'Created dataset from XML element');

for my $attr (qw/ name row_name headers rows attributes /) {
    is_deeply($from_element_ds->$attr, $dataset->$attr, $attr);
}

# aryref
ok(
    !Genome::Report::Dataset->_validate_aryref(
        name => 'data',
        value => undef,
        method => '_validate_aryref',
    ),
    '_validate_aryref failed as expected - no value'
);
ok(
    !Genome::Report::Dataset->_validate_aryref(
        name => 'data',
        value => 'string',
        method => '_validate_aryref',
    ),
    '_validate_string_for_xml failed as expected - not aryref'
);

# xml
ok(
    !Genome::Report::Dataset->_validate_string_for_xml(
        name => 'data',
        value => undef,
        method => '_validate_string_for_xml',
    ),
    '_validate_string_for_xml failed as expected - no value'
);
ok(
    !Genome::Report::Dataset->_validate_aryref(
        name => 'data',
        value => 'string_w_under_scores',
        method => '_validate_string_for_xml',
    ),
    '_validate_string_for_xml failed as expected - value has underscores'
);

is_deeply([$dataset->get_row_values_for_header('assembled')], [5], 'get_row_values_for_header');
ok(!$dataset->get_row_values_for_header(), 'get_row_values_for_header failed as expected - no header');
ok(!$dataset->get_row_values_for_header('not there'), 'get_row_values_for_header failed as expected - header not found');

done_testing();
