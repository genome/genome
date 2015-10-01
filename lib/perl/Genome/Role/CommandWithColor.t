#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More tests => 3;

use Genome::Role::CommandWithColor;

package Genome::Test::CommandWithColor;
use Genome;
class Genome::Test::CommandWithColor {
    is => 'Command::V2',
    roles => 'Genome::Role::CommandWithColor',
};

my $is_in_terminal = 0;
sub _is_running_in_terminal : Overrides(Genome::Role::CommandWithColor) {
    return $is_in_terminal;
}

package main;

subtest basic => sub {
    plan tests => 2;

    my $cmd = Genome::Test::CommandWithColor->create();
    ok($cmd->_status_colors('new'), '_status_colors for known status');
    ok(! $cmd->_status_colors('bogis'), '_status_colors for bogus status');
};

subtest in_terminal => sub {
    _run_terminal_tests($is_in_terminal = 1);
};

subtest not_in_terminal => sub {
    _run_terminal_tests($is_in_terminal = 0);
};

sub _run_terminal_tests {
    my $should_be_colored = shift;

    my $compare = $should_be_colored ? \&like : \&unlike;
    my $test_string = 'hi there';
    my $value = 'the value';
    my %color_map = ( new => 'white' );
    my $ansi_color_regex = qr(\e\[.+m);

    my @tests = (['_color',                 [$test_string, 'white']             => qr($test_string)],
                 ['_colorize_text_by_map',  [$test_string, 'new', %color_map]   => qr($test_string)],
                 ['_color_heading',         [$test_string]                      => qr(===.*$test_string.*===)],
                 ['_color_pair',            [$test_string, $value]              => qr($test_string.*:.*$value)],
                 ['_color_dim',             [$test_string]                      => qr($test_string)],
                );
                
    plan tests => scalar(@tests);
    foreach my $test_data ( @tests ) {
        my($function_name, $args, $output_regex) = @$test_data;
        subtest $function_name => sub {
            plan tests => 2;

            my $cmd = Genome::Test::CommandWithColor->create();
            my $output = $cmd->$function_name(@$args);
            like($output, $output_regex, 'output is formatted properly');
            $compare->($output, $ansi_color_regex, 'produces ansi color if in a terminal');
        };
    }
}

