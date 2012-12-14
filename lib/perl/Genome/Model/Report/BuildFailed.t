#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Test::More;

my $build = Genome::Model::Build->get(107664200); # build for apipe-test-03-MC16s
ok($build, 'Got MC16s build') or die;

my @build_errors = map {
    my $error = Genome::Model::Build::Error->create(
        'build_event_id' => '97901036',
        'stage_event_id' => '97901036',
        'stage' => 'assemble',
        'step' => 'trim-and-screen',
        'step_event_id' => '97901040',
        'error' => 'A really long error message to see if wrapping the text of this error looks good in the report that is generated for users to see why their build failed and what happened to cause it to fail.',
    );
    $error;
} (1..2);
my $generator = Genome::Model::Report::BuildFailed->create(
    build_id => $build->id,
    errors => \@build_errors,
);
ok($generator, 'create');
my $report = $generator->generate_report;
ok($report, 'generate report');

done_testing();
