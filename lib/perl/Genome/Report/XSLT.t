#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Test::More;

use_ok('Genome::Report');
use_ok('Genome::Report::XSLT');

my $dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Report-XSLT';
my $report = Genome::Report->create_report_from_directory($dir.'/Assembly_Stats');
ok($report, 'create report') or die;

#< Valid >#
my $xsl_file = $dir.'/AssemblyStats.txt.xsl';
my $xslt = Genome::Report::XSLT->transform_report(
    report => $report,
    xslt_file => $xsl_file,
);
ok($xslt, 'transformed report');
ok($xslt->{content}, 'Content');
#print $xslt->{content};
ok($xslt->{encoding}, 'Encoding');
is($xslt->{media_type}, 'text/plain', 'Media type');

#< Invalid >#
# no report
my $no_report = Genome::Report::XSLT->transform_report(
    xslt_file => $xsl_file,
);
ok(!$no_report, "Failed as expected, w/o report");
# no xsl_file
my $no_xsl_file = Genome::Report::XSLT->transform_report(
    report => $report,
);
ok(!$no_xsl_file, "Failed as expected, w/o xslt file");

my %media_and_output_types = (
    'application/xml' => 'xml',
    'text/plain' => 'txt',
    'text/html' => 'html',
    rrrr => '',
);

for my $media_type ( keys %media_and_output_types ) {
    is(
        Genome::Report::XSLT->_determine_output_type($media_type), 
        $media_and_output_types{$media_type},
        "Media ($media_type) to output type ($media_and_output_types{$media_type})"
    );
}

done_testing();
exit;

