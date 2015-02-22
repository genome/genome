#!/usr/bin/env genome-perl
use strict;
use warnings;
use above "Genome";
use Test::More;

my $expected_output_dir = $ENV{"GENOME_TEST_INPUTS"} .
  "Genome-Model-ClinSeq-Command-Converge-MutationSpectrum/2014-10-16/";
ok(-e $expected_output_dir, "Found test dir: $expected_output_dir") or die;

my $op_dir = Genome::Sys->create_temp_directory();
ok($op_dir, "Created test dir");

my $op_file = "converge_mutspec.out";
my $bid1="da3b27e756ab49418b23da896698e635";
my $clinseq_build1 = Genome::Model::Build->get($bid1);
ok($clinseq_build1, "Found test build1");
my $bid2="0aeb0be6424446d8a75e42fa377438f2";
my $clinseq_build2 = Genome::Model::Build->get($bid2);
ok($clinseq_build2, "Found test build 2");
my @builds = ($clinseq_build1, $clinseq_build2);
#genome model clin-seq converge mutspec  --builds='model_groups.id=786367aa2edc41e1b4a5d33787a8c003,is_last_complete=1' --outfile=TechD_RNAseq_Metrics.tsv  --verbose=1 --outdir=/tmp/
my $converge_mutspec = Genome::Model::ClinSeq::Command::Converge::MutationSpectrum->create(
      builds => \@builds,
      outfile => $op_file,
      outdir => $op_dir,
      bam_readcount_version => 0.6,
      bq => 0,
      mq => 1,
);
my $return = $converge_mutspec->execute();
ok($clinseq_build2, "Executed succesfully");

my @diff = `diff $expected_output_dir $op_dir`;
ok(@diff == 0, "Found only expected number of differences between expected results and test results")
or do {
  diag("expected: $expected_output_dir\nactual: $op_dir\n");
  diag("differences are:");
  diag(@diff);
  my $diff_line_count = scalar(@diff);
  print "\n\nFound $diff_line_count differing lines\n\n";
  Genome::Sys->shellcmd(cmd => "rm -fr /tmp/last-run-clinseq-converge-mutspec/");
  Genome::Sys->shellcmd(cmd => "mv $op_dir /tmp/last-run-clinseq-converge-mutspec");
};

done_testing();
