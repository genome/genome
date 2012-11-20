#!/usr/bin/env perl
use above "Genome";
use strict;
use warnings;
use Test::More; 

# run with REBUILD on the cmdline to re-generate the reference docs

my @sub_commands = `genome model define 2>&1 | grep '^ ' | awk '{ print \$1 }'`;
for (@sub_commands) { s/[^\w\-]//g; s/\b31m//; s/0m\b//; }
chomp @sub_commands;

my @sub_commands_expected = qw/
  amplicon-assembly
  clin-seq
  convergence
  de-novo-assembly
  differential-expression
  gene-prediction
  genotype-microarray
  germline
  imported-annotation
  imported-assembly
  imported-reference-sequence
  imported-variation-list
  metagenomic-composition-shotgun
  metagenomic-composition16s
  metagenomic-shotgun
  mutational-significance
  phenotype-correlation
  protein-annotation
  reference-alignment
  reference-sequence
  rna-seq
  simple-alignment
  small-rna
  somatic
  somatic-capture
  somatic-validation
  somatic-variation
  test-pipeline
/;

plan tests => ((scalar(@sub_commands_expected)*3)+1);

is("@sub_commands", "@sub_commands_expected", "sub-command list is as expected");

my $expected_dir = __FILE__ . '.expected-output';
my $actual_dir;

if (@ARGV and $ARGV[0] eq 'REBUILD') {
    note("******** regenerating test data in $expected_dir to reset this test case! ***********");
    $actual_dir = $expected_dir;
}
elsif ($ARGV[0]) {
    die "unexpected cmdline options @ARGV: expected nothing or 'REBUILD'";
}
else {
    $actual_dir = Genome::Sys->create_temp_directory;
}

for my $sub_command (@sub_commands) {
    my $actual_out = $actual_dir . '/' . $sub_command;
    my $expected_out = $expected_dir . '/' . $sub_command;

    #my $cmd = "genome model define $sub_command -h >|$actual_out 2>&1";
    #eval { Genome::Sys->shellcmd(cmd => $cmd, allow_failed_exit_code => 1) };
    my $pid = UR::Context::Process->fork();

    if($pid) {
        waitpid($pid,0);
    } else {
        eval {
            local *STDOUT = Genome::Sys->open_file_for_writing($actual_out);
            local *STDERR = *STDOUT;
            local @ARGV = ("model", "define", $sub_command, "-h");
            Genome::Command->execute_with_shell_params_and_exit();
        };

        #the above should exit() so we only get here if it fails somehow
        fail("successful execution of 'genome model define $sub_command -h'");
        exit;
    }

    ok(-s $actual_dir, "output data was generated for $sub_command");

    ok(-s $expected_out, "reference data is available for $sub_command under $expected_dir")
        or next;

    next unless -e $actual_dir;

    my @diff = `diff -r --brief $expected_out $actual_out`;
    is(scalar(@diff), 0, "no differences between actual output and expected output for $sub_command")
        or do {
            for(@diff) { diag($_) };
        }

}

