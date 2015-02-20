#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";
use Test::More;

my $pkg = 'Genome::Model::Tools::Picard::Base';
use_ok($pkg);

class TestPicardCommand {
    is => $pkg,
    has_input => [
        foo => {
            picard_param_name => 'FOO'
        },

        bar => {
            picard_param_name => 'BaR',
            is_many => 1,
        },

        baz => {
            picard_param_name => 'bAz',
            is => 'Boolean',
        },

        qux => {
            picard_param_name => 'QUX',
            is => 'Boolean',
            is_many => 1,
        },
    ]
};


subtest "Command line args" => sub {
    my $obj = TestPicardCommand->create(
        # 1 value
        foo => 'value_foo',
        # 3 more for 4 total
        bar => ['a'..'c'],
        # 5
        baz => 0,
        # we should get 9 params total
        qux => [0, 1, 0, '0 but true'],
        );

    my @args = $obj->_cmdline_args;

    # There are 4 parameters inherited from gmt::picard.
    # we are not interested in testing those right here. they are:
    #   MAX_RECORDS_IN_RAM
    #   VALIDATION_STRINGENCY
    #   TMP_DIR
    #   CREATE_MD5_FILE
    #
    # once refactoring is complete, we might be able to clean this up a little
    is(scalar @args, 9 + 4, 'got the right number of arguments');

    my $extract_values = sub {
        my $prefix = shift;
        return [
            map {my @x = split("=", $_, 2); $x[1]}
            grep {$_ =~ /^$prefix/} @_
            ];
    };

    my %expected = (
        BaR => ['a'..'c'],
        bAz => ['false'],
        FOO => ['value_foo'],
        QUX => ['false', 'true', 'false', 'true'],
        );

    for my $k (keys %expected) {
        my $values = $extract_values->($k, @args);
        is_deeply($values, $expected{$k}, "param $k is correct");
    }
};

done_testing();
