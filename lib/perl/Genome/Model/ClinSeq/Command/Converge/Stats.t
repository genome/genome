#!/usr/bin/env genome-perl
use strict;
use warnings;
use above "Genome";
use Test::More;

my $expected_output_dir = $ENV{"GENOME_TEST_INPUTS"} . "Genome-Model-ClinSeq-Command-Converge-Stats/2014-12-19/";
ok(-e $expected_output_dir, "Found test dir: $expected_output_dir") or die;

my $op_dir = Genome::Sys->create_temp_directory();
ok($op_dir, "Created test dir");

my $op_file = "converge_stats.out";
my $bid1="efbd404356bc4c158ece026c1e018a9b";
my $clinseq_build1 = Genome::Model::Build->get($bid1);
ok($clinseq_build1, "Found test build1");
my $bid2="88bae07356cf4873946b874121772f6d";
my $clinseq_build2 = Genome::Model::Build->get($bid2);
ok($clinseq_build2, "Found test build 2");
my @builds = ($clinseq_build1, $clinseq_build2);
#genome model clin-seq converge stats  --builds='model_groups.id=786367aa2edc41e1b4a5d33787a8c003,is_last_complete=1' --outfile=TechD_RNAseq_Metrics.tsv  --verbose=1 --outdir=/tmp/
my $converge_stats = Genome::Model::ClinSeq::Command::Converge::Stats->create(
      builds => \@builds,
      outdir => $op_dir,
      bam_readcount_version => 0.6,
      bq => 20,
      mq => 20
);
my $return = $converge_stats->execute();
ok($clinseq_build2, "Executed succesfully");

my @diff = `diff $expected_output_dir $op_dir`;
ok(@diff == 0, "Found only expected number of differences between expected results and test results")
or do {
  diag("expected: $expected_output_dir\nactual: $op_dir\n");
  diag("differences are:");
  diag(@diff);
  my $diff_line_count = scalar(@diff);
  print "\n\nFound $diff_line_count differing lines\n\n";
  Genome::Sys->shellcmd(cmd => "rm -fr /tmp/last-run-clinseq-converge-stats/");
  Genome::Sys->shellcmd(cmd => "mv $op_dir /tmp/last-run-clinseq-converge-stats");
};

done_testing();
