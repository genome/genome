#!/usr/bin/env genome-perl

use strict;
use warnings;

package Caller;

use Carp qw(croak);
use Genome::Carp qw(croakf dief);

sub do_croak { croak(@_) }
sub do_croakf { croakf(@_) }

# The two die subs must be on the same line to keep line number the same.
sub do_die { die(@_) }; sub do_dief { dief(@_) }

package main;

use Test::More tests => 2;
use Test::Fatal qw(exception);

my $template = 'Hello %s!';
my @args = ('World');

# The two Caller anonymous subs must be on the same line to keep line number the same.
my ($croak_message, $croakf_message) = map {
    exception { $_->() }
} (sub { Caller::do_croak(sprintf($template, @args)) }, sub { Caller::do_croakf($template, @args) });
is($croakf_message, $croak_message, 'croakf matches croak');

my ($die_message, $dief_message) = map {
    exception { $_->() }
} (sub { Caller::do_die(sprintf($template, @args)) }, sub { Caller::do_dief($template, @args) });
if ($dief_message !~ /\.\n$/) {
    # die/warn end the message with a period, but earlier versions of Carp::croak did not.
    chomp $dief_message;
    $dief_message .= ".\n";
}
is($dief_message, $die_message, 'dief correctly approximated die');
