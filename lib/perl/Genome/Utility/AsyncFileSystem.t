#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";
use Test::More;

plan tests => 28;

use_ok 'Genome::Utility::AsyncFileSystem';

use Genome::Utility::AsyncFileSystem qw(on_each_line);

# on_each_line

on_each_line_checker(
    'simple lines',
    [ "foo\n", "bar\n", "baz\n", "boo\n" ],
    [ "foo\n", "bar\n", "baz\n", "boo\n" ]
);

on_each_line_checker(
    'partial lines',
    [ "foo",   "\nbar\nbaz", "\nboo", "\n" ],
    [ "foo\n", "bar\n",      "baz\n", "boo\n" ]
);

on_each_line_checker(
    'without trailing \n',
    [ "foo", "\nbar\nbaz", "\nboo" ],
    [ "foo\n", "bar\n", "baz\n", "boo" ]
);

on_each_line_checker( 'without any \n',
    ["foo bar baz boo"], ["foo bar baz boo"] );

on_each_line_checker( 'input chunk', ["foo\nbar\nbaz\nboo"],
    [ "foo\n", "bar\n", "baz\n", "boo" ] );

sub on_each_line_checker {
    my ( $prefix, $inputs, $expected ) = @_;
    my $line_closure = on_each_line {
        is( $_[0], shift @$expected, $prefix . ' on_each_line called arg' );
    };

    for (@$inputs) {
        $line_closure->($_);
    }
    $line_closure->();

    if (@$expected) {
        fail( $prefix . ' expected array exhausted' );
    } else {
        pass( $prefix . ' expected array exhausted' );
    }
}

