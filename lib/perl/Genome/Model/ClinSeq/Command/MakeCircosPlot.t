#!/usr/bin/env genome-perl
use strict;
use warnings;
use above "Genome"; # only use "above" in test cases, never in modules
use Genome::Model::ClinSeq::TestData;
use Test::More tests => 5;

# this gets a canonical test build, and is set up to not really hit the db
# we temporarily have it overridden to test against Obi's example
my $base_dir= $ENV{"GENOME_TEST_INPUTS"} . "Genome-Model-ClinSeq-Command-MakeCircosPlot/2013-10-22";
my $expected_output_dir = "$base_dir/expected-output";
ok(-e $expected_output_dir, "expected output dir exists: $expected_output_dir");

#my $test_ids = Genome::Model::ClinSeq::TestData::load(base_dir => "$base_dir/input_dir");
#my $test_build_id = $test_ids->{CLINSEQ_BUILD};
my $test_build_id = "aee9a31051754702a9c2835d63abc812"; ## temp override
my $test_build = Genome::Model::Build->get($test_build_id);
ok($test_build, "got test build $test_build_id");

# this directory lives under /gsc/var/cache/testsuite/data at TGI
#/gscmnt/gc13003/info/test_suite_data//Genome-Model-ClinSeq-Command-MakeCircosPlot/2013-10-01/expected-output

# make a temp directory for output
my $actual_output_dir = Genome::Sys->create_temp_directory();
ok(-e $actual_output_dir, "test output dir exists");

# run the command
my $result = Genome::Model::ClinSeq::Command::MakeCircosPlot->execute(
    build => $test_build,
    output_directory => $actual_output_dir,
);
print $result."\n";
ok($result, "execution succeeded");

# verify results


    my @differences = `diff $expected_output_dir $actual_output_dir`;
    is(scalar(@differences), 35, "only expected differences found: diff $expected_output_dir $actual_output_dir")
        or do {
            #print "DIFF:\n", @differences;
            # un-comment these to keep the output for debugging 
            my $debug_location = "/tmp/last-output-make-circos-plot-$$";
            system "mv $actual_output_dir $debug_location";
        };

#TODO create tests for WGS only, WGS+Exome, and WGS+RNA. now this only tests for all three
