#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";

require File::Compare;
use Test::More;

use_ok( 'Genome::Model::Tools::Velvet::Metrics' ) or die;

my $data_dir = $ENV{GENOME_TEST_INPUTS} . "/Genome-Model-Tools-Velvet/Metrics/v8";
my $example_metrics = $data_dir.'/metrics.txt';

#create temp test dir
my $temp_dir = Genome::Sys->create_temp_directory();
my $metrics_file = $temp_dir.'/metrics.txt';

#create stats
my $metrics = Genome::Model::Tools::Velvet::Metrics->create(
    assembly_directory => $data_dir,
    output_file => $metrics_file,
    min_contig_length => 1,
);
ok($metrics, "create");
$metrics->dump_status_messages(1);
ok($metrics->execute, "execute");

is(File::Compare::compare($metrics_file, $example_metrics), 0, 'metrics file matches');

#print "gvimdiff $metrics_file $example_metrics\n"; system "gvimdiff $metrics_file $example_metrics"; <STDIN>;
done_testing();
exit;

