#!/usr/bin/env genome-perl

use strict;
use warnings FATAL => qw(all);

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
};

use above 'Genome';
use Test::More tests => 6;

my $TEST_STATUS = 'unique-test-status';

{
    package Genome::Config::AnalysisProject::Command::Test;
    class Genome::Config::AnalysisProject::Command::Test {
        is => 'Genome::Config::AnalysisProject::Command::Base',
    };

    sub valid_statuses {
        return ($TEST_STATUS);
    }
}

my $test_class = 'Genome::Config::AnalysisProject::Command::Test';


my $help_text = $test_class->help_usage_complete_text;
ok($help_text, 'can generate help text');
like($help_text, qr/\Q$TEST_STATUS\E/, 'found test status in documentation');

my $anp = Genome::Config::AnalysisProject->__define__(
    name => 'test anp for base command test',
    status => 'Hold',
);

my $cmd = $test_class->create(analysis_project => $anp);
isa_ok($cmd, 'Genome::Config::AnalysisProject::Command::Base', 'test command');

my @errors = $cmd->__errors__;
is(scalar(@errors), 1, 'got error about invalid status');
is($errors[0]->desc, qq(Can't test using analysis project with status: Hold), 'error text is identically what we thought it would be');

$anp->status($TEST_STATUS);

@errors = $cmd->__errors__;
is(scalar(@errors), 0, 'no error with expected status');
