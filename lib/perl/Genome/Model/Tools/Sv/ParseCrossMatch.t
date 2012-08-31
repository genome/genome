#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";
use Test::More;
use File::Compare;
use Data::Dumper;

plan tests => 7;

use_ok( 'Genome::Model::Tools::Sv::ParseCrossMatch');

my $test_input_dir  = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Sv-ParseCrossMatch/';
my $cm_file = $test_input_dir . 'cm.file';

my $cm = Genome::Model::Tools::Sv::ParseCrossMatch->create(
    input_file => $cm_file,
);

ok($cm, 'ParseCrossMatch created ok');

my $pos = $cm->dcpos->{458}->[0];

is($pos->{score}, '354', 'score matches');
is($pos->{type}, 'S', 'type matches');
is($pos->{base}, 'C', 'base matches');
is($pos->{rpos}, '194', 'rpos matches');
is($pos->{aln_orient}, 'C', 'aln_orient matches');


