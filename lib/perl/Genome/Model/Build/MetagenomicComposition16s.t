#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Genome::Model::Build::MetagenomicComposition16s::TestBuildFactory;
use Test::More;

use_ok('Genome::Model::Build::MetagenomicComposition16s') or die;

diag('454');
my ($build, $example_build) = Genome::Model::Build::MetagenomicComposition16s::TestBuildFactory->build_with_example_build_for_454;
ok($build, 'Got mc16s 454 build');
is($build->calculate_estimated_kb_usage, 2000, 'Estimated kb usage');
my @amplicon_sets = $build->amplicon_sets;
is(@amplicon_sets, 3, 'Amplicon sets');
@amplicon_sets = $build->amplicon_sets_for_processing;
is(@amplicon_sets, 4, 'Amplicon sets for processing');

diag('Sanger');
($build, $example_build) = Genome::Model::Build::MetagenomicComposition16s::TestBuildFactory->build_with_example_build_for_sanger;
ok($build, 'Got mc16s sanger build');
is($build->calculate_estimated_kb_usage, 30_000, 'kb usage');
@amplicon_sets = $build->amplicon_sets;
is(@amplicon_sets, 1, 'Amplicon sets');
@amplicon_sets = $build->amplicon_sets_for_processing;
is(@amplicon_sets, 1, 'Amplicon sets for processing');

done_testing();
