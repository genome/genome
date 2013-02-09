#!/usr/bin/env perl
use strict;
use warnings;
use above "Genome";
use Test::More tests => 8;
use Genome::Model::Tools::Htseq::Count;

$ENV{UR_DBI_NO_COMMIT} = 1;

my $test_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Htseq-Count';

# find the test alignment result

#my $b = Genome::Model::Build->get(133351985);
#ok($b, "got test build " . $b->__display_name__);

my $a = Genome::InstrumentData::AlignmentResult->get(133352072);
ok($a, "got alignment result " . $a->__display_name__);

# The tool works with the AR output_dir by default
# but will take a shortcut and use the scratch sorted BAM if it exists.
# This means the tool can be run after the AR is complete,
# or while it is being generated, and in the later case will run more efficiently.

# Swith the temp scratch directory on the AR to the one we have set up
# it only has chromosome 22 data so it will be faster.  This will make the tool
# believe it is running during alignment and we will test its more efficient execution.

$a->temp_scratch_directory($test_dir . '/fake-scratch-dir-for-alignment-result');
ok(-d $a->temp_scratch_directory, "set the temp_scratch_directory for the software result to our test data to " . $a->temp_scratch_directory);


# verify test inputs and outputs exist

my $input_bam = $a->temp_scratch_directory . '/raw_all_sequences.bam.sort.bam';
ok(-e $input_bam, "found input sorted test bam which only has chr22");

my $expected_out1 = $test_dir . '/expected-outputs/gene-counts.tsv';
ok(-e $expected_out1, "found comparison output file 1: $expected_out1");

my $expected_out2 = $test_dir . '/expected-outputs/transcript-counts.tsv'; 
ok(-e $expected_out2, "found comparison output file 2: $expected_out2");

# run in /tmp
#my $test_outdir = Genome::Sys->create_temp_directory();
#ok(-d $test_outdir, "created test output directory");

# before running, ensure results do not exist previously
my $test_name = $ENV{GENOME_SOFTWARE_RESULT_TEST_NAME} ||= "testsuite " . UR::Context->now . " " . Sys::Hostname::hostname() . "-$$.";
my $result_exists = Genome::Model::Tools::Htseq::Result->get(
    alignment_result => $a,
    test_name => $test_name
);
ok(!$result_exists, "no result already in the system for this test") or die "contact informatics!";

my $new_result = Genome::Model::Tools::Htseq::Count->execute(
    alignment_results => [$a],
    #output_dir => $test_outdir, # remove when automatic SR generateion is in place
    app_version => '0.5.3p9',
    results_version => 1,
    limit => 2000,
);

# diff results

for my $pair(
    [$new_result->output_dir . '/gene-counts.tsv',          $test_dir . '/expected-outputs/gene-counts.tsv'],
    [$new_result->output_dir . '/transcript-counts.tsv',    $test_dir . '/expected-outputs/transcript-counts.tsv'],
) {
    my ($actual,$expected) = @$pair;
    my @diff = `diff $actual $expected`;
    ok(scalar(@diff)==0, "no differences for $actual vs $expected")
        or do {
            if ($ENV{USER} ne 'apipe-tester') {
                my $dir = '/tmp/last-htseq-failed-test-' . $$;
                diag("moving failed results to $dir");
                rename $new_result->output_dir, $dir;
            }
        };
}

my $found_result_after = Genome::Model::Tools::Htseq::Result->get(
    alignment_result => $a,
    test_name => $test_name
);
ok($found_result_after, "found a result after running the command");
is($found_result_after, $new_result, "it matches the result returned from the tool");

# remove the new result to make sure removal works
# because we have no-commit turned on it will be removed at test exit regardless.
$found_result_after->delete;
ok($found_result_after->isa("UR::DeletedRef"), "deletion worked");


