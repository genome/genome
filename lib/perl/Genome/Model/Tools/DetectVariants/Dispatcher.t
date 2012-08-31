#!/usr/bin/env genome-perl

use strict;
use warnings;

use File::Path;
use File::Temp;
use Test::More tests => 7;
#use Test::More skip_all => 'test in development';
use above 'Genome';

#Use this check if tests are added to actually run detectors
#my $archos = `uname -a`;
#if ($archos !~ /64/) {
#    plan skip_all => "Must run from 64-bit machine";
#} else {
#    plan tests => 5;
#}

#Parsing tests
my $detector_string = '(samtools && varscan) || maq';

my $dispatcher = 'Genome::Model::Tools::DetectVariants::Dispatcher';

my $tree = $dispatcher->parse_detector_string($detector_string);
ok($tree, 'able to parse detector string');


#Test out the tree walker by using it to reconstruct a description of detectors to run from the detector string
my $branch_case = sub {
    my $self = shift;
    my ($combination, $subtrees, $branch_case, $leaf_case) = @_;
    
    return '(' . join(' ', $self->walk_tree($subtrees->[0], $branch_case, $leaf_case), $combination, $self->walk_tree($subtrees->[1], $branch_case, $leaf_case)) . ')';
};

my $leaf_case = sub {
    my $self = shift;
    my ($detector_name, $index, $branch_case, $leaf_case) = @_;
    
    return $detector_name . ' (call #' . $index . ')';
};

my $parsed_string = $dispatcher->walk_tree($tree, $branch_case, $leaf_case);

my $expected_parsed_string = '((samtools (call #0) intersect varscan (call #1)) union maq (call #2))';
is($parsed_string, $expected_parsed_string, 'tree-walker returns expected result');

my $varscan_class = $dispatcher->detector_class('varscan');
is($varscan_class, 'Genome::Model::Tools::DetectVariants::Varscan', 'resolved detector class for varscan');
ok($dispatcher->is_valid_detector($varscan_class), 'determined that varscan is in fact a detector');

my $nonexistent_class = $dispatcher->detector_class('non-existent');
ok((not $dispatcher->is_valid_detector($nonexistent_class)), 'determined that non-existent is not a detector');

my $somatic_detector_string = '(somatic sniper || somatic varscan)';
my $somatic_tree = $dispatcher->parse_detector_string($somatic_detector_string);
ok($somatic_tree, 'able to parse detector_string');

my $somatic_parsed_string = $dispatcher->walk_tree($somatic_tree, $branch_case, $leaf_case);

my $expected_somatic_parsed_string = '(somatic sniper (call #0) union somatic varscan (call #1))';
is($somatic_parsed_string, $expected_somatic_parsed_string, 'tree-walker returns expected somatic result');
