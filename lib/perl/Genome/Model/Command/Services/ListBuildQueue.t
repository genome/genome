#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::More tests => 4;

use Genome;

use Genome::Model::Command::Services::ListBuildQueue;
use Genome::Test::Factory::Model::ReferenceAlignment;

my @m = map { Genome::Test::Factory::Model::ReferenceAlignment->setup_object() } 1..3;
for (@m) { $_->build_requested(1) }
my @model_ids = map { $_->id } @m;
do {
    no warnings qw(once redefine);
    *Genome::Model::create_iterator = sub {
        package Genome::Model;
        my $self = shift;
        return $self->SUPER::create_iterator(id => \@model_ids, @_);
    };
};

do {
    my $line_count = 0;
    no warnings qw(once redefine);
    local *Genome::Model::Command::Services::ListBuildQueue::print_model_ids = sub { $line_count++ };

    Genome::Model::Command::Services::ListBuildQueue->execute(
        count => 1,
    );
    is($line_count, scalar(@m));
};

do {
    my $line_count = 0;
    no warnings qw(once redefine);
    local *Genome::Model::Command::Services::ListBuildQueue::print_model_ids = sub { $line_count++ };

    Genome::Model::Command::Services::ListBuildQueue->execute(
        count => 2,
    );
    my $expected = int(@m / 2) + @m % 2;
    is($line_count, $expected);
};

do {
    my $line_count = 0;
    no warnings qw(once redefine);
    my $scheduled_builds_for = 0;
    local *Genome::Model::Command::Services::ListBuildQueue::scheduled_builds_for = sub { $scheduled_builds_for++ };
    local *Genome::Model::Command::Services::ListBuildQueue::print_model_ids = sub { $line_count++ };

    Genome::Model::Command::Services::ListBuildQueue->execute(
        max => 1,
    );
    is($line_count, 1);
};

do {
    my $output;
    open my $h, '>', \$output;
    my $f = \&Genome::Model::Command::Services::ListBuildQueue::print_model_ids;
    $f->($h, "\0", '123');
    chomp $output;
    is($output, "123\0", "must be null terminated arguments even for one argument");
};
