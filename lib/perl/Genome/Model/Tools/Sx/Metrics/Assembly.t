#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";

require File::Compare;
use Test::More;

use_ok('Genome::Model::Tools::Sx::Metrics::Assembly') or die;

#check testsuite data files
my $data_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Sx/Metrics/Assembly/v6';
ok(-d $data_dir, "Data dir exists") or die;

#create temp test dir
my $temp_dir = Genome::Sys->create_temp_directory();
ok(-d $temp_dir, "Created temp test dir");

my $metrics = Genome::Model::Tools::Sx::Metrics::Assembly->create(
    major_contig_threshold => 300,
    tier_one => 3550,
    tier_two => 3550,
);
ok($metrics, "create") or die;
$metrics->add_contigs_file($data_dir.'/contigs.bases:type=fasta');
for my $reads_file ( $data_dir.'/1_fastq', $data_dir.'/2_fastq' ) {
    $metrics->add_reads_file($reads_file.':type=sanger');
}

#print $metrics->to_xml;
my $text = $metrics->transform_xml_to('txt');
ok($text, 'got text');
my $metrics_file = $temp_dir.'/metrics.txt';
my $fh = Genome::Sys->open_file_for_writing($metrics_file);
$fh->print("$text");
$fh->close;

my $example_metrics = $data_dir.'/metrics.txt';
is(File::Compare::compare($metrics_file, $example_metrics), 0, "files match");

#print "gvimdiff $metrics_file $example_metrics\n"; system "gvimdiff $metrics_file $example_metrics\n"; <STDIN>;
done_testing();
