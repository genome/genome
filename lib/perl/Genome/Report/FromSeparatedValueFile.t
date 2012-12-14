#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Test::More;

use_ok('Genome::Report::FromSeparatedValueFile') or die;

my $svr = Genome::Utility::IO::SeparatedValueReader->create(
    input => $ENV{GENOME_TEST_INPUTS} . '/Genome-Utility-IO/albums.csv',
);
ok($svr, 'create svr') or die;

my %params = (
    name => 'Report from Albums SVF',
    description => 'Albums on Hand Today',
    svr => $svr,
);

my $generator = Genome::Report::FromSeparatedValueFile->create(%params);
ok($generator, 'create generator');
can_ok($generator, '_add_to_report_xml');

my $report = $generator->generate_report;
ok($report, 'Generated report');
ok($report->xml_string, 'report xml string');

for my $attr (qw/ name description svr /) {
    my $value = delete $params{$attr};
    ok(!Genome::Report::FromSeparatedValueFile->create(%params), "failed to create w/o $attr");
    $params{$attr} = $value;
}

done_testing();
