#!/usr/bin/env genome-perl

use strict;
use warnings;
    
use above 'Genome';

use Test::More;

use_ok('Genome::Report::Email');

my $test_report_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Report/Build_Start';
    
my $report = Genome::Report->create_report_from_directory($test_report_dir);
ok($report, 'create report') or die;


my %valid_params = (
    report => $report,
    to => [Genome::Config->user_email], # can be string or aryref
    xsl_files => [ $report->generator->get_xsl_file_for_html ],
);

#< Valid >#
my $valid = Genome::Report::Email->send_report(%valid_params);
ok($valid, 'Sent report');

#< Invalid >#
for my $attr (qw/ to report xsl_files /) {
    my $val = delete $valid_params{$attr};
    my $invalid = Genome::Report::Email->send_report(%valid_params);
    ok(!$invalid, 'Failed as expected - no '.$attr);
    $valid_params{$attr} = $val;
}

done_testing();
