#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::More;
use File::Temp;

use above "Genome";

test_stderr_tie();

done_testing();

sub test_stderr_tie {
    my $tmp_file = File::Temp->new();

    UR::Object::Type->define(
        class_name => 'TestLogger',
        roles => 'Genome::Role::Logger',
    );
    my $logger = TestLogger->create(
        screen => 0,
        log_file => "$tmp_file",
        tie_stderr => 1,
    );
    $logger->delegate_logger();

    my $test_message = "Test.\n";
    my $expected_file_content = "STDERR: $test_message";
    print STDERR $test_message;

    my $tmp_file_content = <$tmp_file>;
    is($tmp_file_content, $expected_file_content, 'tmp_file content matches expected_file_content');
}
