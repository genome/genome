#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";
use Genome::Test::Factory::AnalysisProject;

use Test::More tests => 6;

use_ok('Genome::Config::Tag::AnalysisProject::SubjectMapping');

my $tag = Genome::Config::Tag->create(
    name => 'testing_linkage_to_subject_mappings',
    description => 'this is only a test',
);
isa_ok($tag, 'Genome::Config::Tag');

my $test_anp = Genome::Test::Factory::AnalysisProject->setup_object();
my $subject_mapping = Genome::Config::AnalysisProject::SubjectMapping->create(
    analysis_project => $test_anp,
);
isa_ok($subject_mapping, 'Genome::Config::AnalysisProject::SubjectMapping');

my $tag_assignment = Genome::Config::Tag::AnalysisProject::SubjectMapping->create(
    tag => $tag,
    subject_mapping => $subject_mapping,
);
isa_ok($tag_assignment, 'Genome::Config::Tag::AnalysisProject::SubjectMapping');

is($tag->subject_mappings, $subject_mapping, 'tag sees new assignment');
is($subject_mapping->tags, $tag, 'subject mapping sees new assignment');
