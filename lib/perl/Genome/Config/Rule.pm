package Genome::Config::Rule;

use strict;
use warnings;

use Genome;

class Genome::Config::Rule {
    is => 'UR::Object',
    is_transactional => 0,
    doc => 'encapsulates a single rule',
    has => [
        method_chain => {
            is => 'ARRAY'
        },
        expected_value => {
            is => 'String'
        }
    ]
};

sub create_from_hash {
    my $class = shift;
    my $rules_hash = shift;

    die('Need to pass in a hash ref of rules!')
        unless $rules_hash && ref $rules_hash eq 'HASH';

    my @rule_objects;
    while (my ($key, $val) = each %$rules_hash) {
        my $method_chain = $class->_parse_method_chain($key);
        my $expected_value = $val;
        push @rule_objects, $class->create(
            method_chain => $method_chain,
            expected_value => $expected_value
        );
    }

    return @rule_objects;
}

sub _parse_method_chain {
    my $class = shift;
    my $method_string = shift;

    die('No method string given!') unless $method_string;

    return [split('->', $method_string)];
}

1;
