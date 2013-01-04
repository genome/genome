#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More tests => 5;

use Genome::Sys;

my $tmp = Genome::Sys->create_temp_directory();

# TODO: use a "testing" cDNA model to get build rather than this real build
my $build_id = 98314469;

my $r = Genome::Model::ReferenceAlignment::Report::ReferenceCoverage->create(
    build_id => $build_id,
);
ok($r, "created a new report");

my $v = $r->generate_report;
ok($v, "generation worked");

my $result = $v->save($tmp);
ok($result, "saved to $tmp");

my $name = $r->name;
$name =~ s/ /_/g;

ok(-d "$tmp/$name", "report directory $tmp/$name is present");
ok(-e "$tmp/$name/report.xml", 'xml report is present');
