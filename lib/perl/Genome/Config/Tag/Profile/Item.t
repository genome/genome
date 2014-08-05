#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";
use Genome::Test::Factory::AnalysisProject;

use Test::More tests => 4;

use_ok('Genome::Config::Tag::Profile::Item');

my $tag = Genome::Config::Tag->create(
    name => 'testing_linkage_to_profile_items',
    description => 'this is only a test',
);
isa_ok($tag, 'Genome::Config::Tag');

my $test_anp = Genome::Test::Factory::AnalysisProject->setup_object();
my $profile_item = Genome::Config::Profile::Item->create(
    analysis_project => $test_anp,
);
isa_ok($profile_item, 'Genome::Config::Profile::Item');

my $tag_assignment = Genome::Config::Tag::Profile::Item->create(
    tag => $tag,
    profile_item => $profile_item,
);
isa_ok($tag_assignment, 'Genome::Config::Tag::Profile::Item');

