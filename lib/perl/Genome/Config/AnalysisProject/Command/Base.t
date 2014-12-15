#!/usr/bin/env genome-perl

use strict;
use warnings FATAL => qw(all);

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
};

use above 'Genome';
use Test::More tests => 3;

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

isa_ok($test_class, 'Genome::Config::AnalysisProject::Command::Base', 'test command');

my $help_text = $test_class->help_usage_complete_text;
ok($help_text, 'can generate help text');
like($help_text, qr/\Q$TEST_STATUS\E/, 'found test status in documentation');

