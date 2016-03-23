#!/usr/bin/env genome-perl
use strict;
use warnings;
use above "Genome";
use Test::More;

my $expected_output_dir =
    Genome::Config::get('test_inputs') . "Genome-Model-ClinSeq-Command-Converge-CancerRelevanceScore/2014-09-04/";
ok(-e $expected_output_dir, "Found test dir: $expected_output_dir") or die;

my $op_dir = Genome::Sys->create_temp_directory();
ok($op_dir, "Created test dir");

my $op_file        = "converge_cancerRelevanceScore.out";
my $bid1           = "d0f2169786074f6e9609a9f6de3754c7";
my $clinseq_build1 = Genome::Model::Build->get($bid1);
ok($clinseq_build1, "Found test build1");
my $bid2           = "638b5e3c3e7548df94351907d6d8c714";
my $clinseq_build2 = Genome::Model::Build->get($bid2);
ok($clinseq_build2, "Found test build 2");
my @builds = ($clinseq_build1, $clinseq_build2);

#genome model clin-seq converge cancer-relevance-score
#   --builds='model_groups.id=786367aa2edc41e1b4a5d33787a8c003,is_last_complete=1'
#   --outdir=/tmp/
my $converge_cancerRelevanceScore = Genome::Model::ClinSeq::Command::Converge::CancerRelevanceScore->create(
    builds                => \@builds,
    outfile               => $op_file,
    outdir                => $op_dir,
    bam_readcount_version => 0.6,
);

my $return = $converge_cancerRelevanceScore->execute();
ok($clinseq_build2, "Executed succesfully");

my @diff = `diff $expected_output_dir $op_dir`;
my $ok = ok(@diff == 0, "Found only expected number of differences between expected results and test results");
unless ($ok) {
    diag("expected: $expected_output_dir\nactual: $op_dir\n");
    diag("differences are:");
    diag(@diff);
    my $diff_line_count = scalar(@diff);
    print "\n\nFound $diff_line_count differing lines\n\n";
    Genome::Sys->shellcmd(cmd => "rm -fr /tmp/last-run-clinseq-converge-cancerRelevanceScore/");
    Genome::Sys->shellcmd(cmd => "mv $op_dir /tmp/last-run-clinseq-converge-cancerRelevanceScore");
}

done_testing();
