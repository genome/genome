#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";

use Data::Dumper 'Dumper';
use Test::More;
use XML::LibXML;

use_ok('Genome::Report');
my $dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Report/';
my $reports_dir = $dir.'/xml_reports';
my $tmp_dir = Genome::Sys->base_temp_directory;

my $xml = XML::LibXML->new->parse_string(<<EOS
<?xml version="1.0"?>
<report>
  <datasets>
    <stats>
      <stat>
        <assembled>4</assembled>
        <attempted>5</attempted>
      </stat>
    </stats>
  </datasets>
   <report-meta>
    <name>Assembly Stats</name>
    <description>Assembly Stats for Amplicon Assembly (Name &lt;mr. mock&gt; Build Id &lt;-10000&gt;)</description>
    <date>2009-05-29 10:19:10</date>
    <generator>Genome::Model::AmpliconAssembly::Report::AssemblyStats</generator>
    <generator-params>
      <build-id>-10000</build-id>
      <amplicons>HMPB-aad16a01</amplicons>
      <amplicons>HMPB-aad16c10</amplicons>
    </generator-params>
  </report-meta>
</report>
EOS
);
ok($xml, 'parse xml string');
my $report = Genome::Report->create(xml => $xml);
ok($report, 'create report');

# meta
my %report_meta = (
    name => 'Assembly Stats',
    description => 'Assembly Stats for Amplicon Assembly (Name <mr. mock> Build Id <-10000>)',
    date => '2009-05-29 10:19:10',
    generator => 'Genome::Model::AmpliconAssembly::Report::AssemblyStats',
    generator_params => {
        build_id => [ -10000 ],
        amplicons => [qw/ HMPB-aad16a01 HMPB-aad16c10 /],
    },
);
for my $attr ( keys %report_meta ) {
    is_deeply($report->$attr, $report_meta{$attr}, $attr);
}

# dir
ok(!$report->parent_directory, 'parent directory is undef');
ok(!$report->directory, "Can't access directory w/o setting parent directory");

# datasets
my @datasets = $report->get_dataset_nodes;
ok(@datasets, 'datasets') or die;
my ($stats) = $report->get_dataset_nodes_for_name('stats');
ok($stats, 'Dataset node: stats') or die;
is($datasets[0]->nodeName, $stats->nodeName, 'Dataset name: stats');

# Save - fails
ok(!$report->save, 'Failed as expected - no directory');
ok(!$report->save('invalid_directory'), 'Failed as expected - invalid directory');

# Save
ok($report->save( $tmp_dir ), 'Saved report to tmp dir');
is($report->directory, $tmp_dir.'/Assembly_Stats', 'Directory name matches');

# Resave 
ok(!$report->save( $tmp_dir ), 'Failed as expected - resave report');
ok($report->save($tmp_dir, 1), 'Overwrite report');

# Get it w/ parent dir
my @reports = Genome::Report->create_reports_from_parent_directory( $reports_dir );
is(@reports, 2, 'Got two reports using create_reports_from_parent_directory');
is($reports[0]->name, 'Test Report 1', 'Report 1 has correct name');
is($reports[0]->parent_directory, $reports_dir, 'Report 1 has correct parent directory');
is($reports[1]->name, 'Test Report 2', 'Report 2 has correct name');
is_deeply(
    $reports[1]->generator_params,
    { 'scalar' => [qw/ yes /], array => [qw/ 1 2 /] },
    'Test Report 1 - gen params'
);

my $valid_name = 'Test Report';
my $valid_dir = $dir;
my @reports;

# get - can't
eval {
    @reports = Genome::Report->get(xml => _xml_string());
};
print "$@\n";
ok(!@reports, 'Failed as expected - get');

# create w/ invalid parent dir
eval {
    @reports = Genome::Report->create_reports_from_parent_directory('no_way_this_dir_exists');
};
print "$@\n";
ok(!@reports, 'Failed as expected - create reports w/ invalid parent_directory');

done_testing();
