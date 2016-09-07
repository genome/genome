#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::More tests => 2;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
};

use above 'Genome';

use Genome::Model::Command::Services::ListBuildQueue;
use Genome::Test::Factory::Model::ReferenceAlignment;
use Genome::Test::Factory::AnalysisProject;

use List::MoreUtils qw(uniq);

my @m = setup();

do {
    my $max = scalar(@m) - 1;
    my $count;
    no warnings qw(once redefine);
    local *Genome::Model::Command::Services::ListBuildQueue::run_as_list = sub { uniq map { $_->run_as } @m };
    local *Genome::Model::Command::Services::ListBuildQueue::builds_for = sub { 0 };
    local *Genome::Model::Command::Services::ListBuildQueue::print = sub { $count++ };
    Genome::Model::Command::Services::ListBuildQueue->execute(
        max => $max,
    );
    is($count, $max, 'max limited output');
};

do {
    my $c = Genome::Model::Command::Services::ListBuildQueue->create();
    my $got = $c->sprint('123');
    is($got, '123' . $c->delimiter, 'delimiter is appended');
};

sub setup {
    my $anp = Genome::Test::Factory::AnalysisProject->setup_object();
    my @m = map { Genome::Test::Factory::Model::ReferenceAlignment->setup_object() } 1..3;
    for (@m) {
        $anp->add_model_bridge(model_id => $_->id);
        $_->build_requested(1);
    }
    my @model_ids = map { $_->id } @m;
    do {
        no warnings qw(once redefine);
        *Genome::Model::create_iterator = sub {
            package Genome::Model;
            my $self = shift;
            return $self->SUPER::create_iterator(id => \@model_ids, @_);
        };
    };
    return @m;
}
