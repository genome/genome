#!/usr/bin/env perl
use strict;
use warnings;
use above "Genome";
use Test::More tests => 9;
use Genome::Model::Tools::Htseq::Count;

$ENV{UR_DBI_NO_COMMIT} = 1;

# try with a pair and verify that merging results works
# gmt htseq count --align id:133654083/133965883 --app-version 0.5.4p1

my $test_name = $ENV{GENOME_SOFTWARE_RESULT_TEST_NAME} ||= "testsuite " . UR::Context->now . " " . Sys::Hostname::hostname() . "-$$.";

my @alignment_result_ids = (133654083,133965883);
my @alignment_results = Genome::InstrumentData::AlignmentResult->get(id => \@alignment_result_ids);
my $set = Genome::InstrumentData::AlignmentResult->define_set(id => \@alignment_result_ids);

my @sr = grep { ! $_->test_name } Genome::Model::Tools::Htseq::Count::Result->get("alignment_results.id in" => \@alignment_result_ids);
is(scalar(@sr), 3, "got three existing results") or die; 

# put our test name here temporarily so we don't regenerate the underlying results
for (@sr) {
    $_->set_test_name($test_name);
}

# generate a composite result
my $cmd = Genome::Model::Tools::Htseq::Count->execute(
    alignment_results => \@alignment_results,
    result_version => 1,
);
ok($cmd, "executed command with two inputs");
my $result = $cmd->result();
ok($result, "got result");

ok(-d $result->output_dir . '/underlying_results', "found dir of underlying results");

my @links = glob($result->output_dir() . '/underlying_results/*');
is(scalar(@links), 2, "found two links under it")
    or diag(join("\n",@links));

my @link_dirs = map { Cwd::abs_path($_) } @links;
note "link dirs is @link_dirs";
my %link_dirs = map { $_ => 1 } @link_dirs;

my @use_associations = Genome::SoftwareResult::User->get(user_id => $result->id);
is(scalar(@use_associations),2, "linked to two other results");

my @uses = map { $_->software_result } @use_associations;
is(scalar(@uses), 2, "uses to other results");

for my $used_sr (@uses) {
    my $dir = $used_sr->output_dir;
    ok($link_dirs{$dir}, "found output dir $dir in links @link_dirs");
}


