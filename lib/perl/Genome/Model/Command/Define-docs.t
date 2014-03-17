#!/usr/bin/env genome-perl
use above "Genome";
use strict;
use warnings;
use Test::More;
use Command::Shell;

Genome::SoftwareResult->class; #This emits warnings under no-commit, so get them out of the way.


# run with REBUILD on the cmdline to re-generate the reference docs
# run with REBUILD $typename on the cmdline to re-generate for just one pipeline

my @sub_commands = `genome model define 2>&1 | grep '^ ' | awk '{ print \$1 }'`;
for (@sub_commands) { s/[^\w\-]//g; s/\b31m//; s/0m\b//; }
chomp @sub_commands;

my @sub_commands_expected = qw/
  amplicon-assembly
  clin-seq
  comparison
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

plan tests => ((scalar(@sub_commands_expected)*5)+1);

is("@sub_commands", "@sub_commands_expected", "sub-command list is as expected");

my $expected_dir = __FILE__ . '.expected-output';
my $actual_dir;
my $rebuild = 0;

if (@ARGV and $ARGV[0] eq 'REBUILD') {
    note("******** regenerating test data in $expected_dir to reset this test case! ***********");
    $actual_dir = $expected_dir;
    $rebuild = 1;
    if ($ARGV[1]) {
        shift @ARGV;
        @sub_commands = @ARGV;
        note("******* regenerationg just for pipelines: @sub_commands ***********");
    }
}
elsif ($ARGV[0]) {
    die "unexpected cmdline options @ARGV: expected nothing or 'REBUILD', got " . $ARGV[0];
}
else {
    $actual_dir = Genome::Sys->create_temp_directory;
}

for my $sub_command (@sub_commands) {
    my $module = module_from_sub_command($sub_command);
    use_ok($module, "Can use $module");

    my $actual_out = $actual_dir . '/' . $sub_command;
    my $expected_out = $expected_dir . '/' . $sub_command;

    eval {
        if ($rebuild) {
            # the previous results may already exist, which will fail the open below
            note("Removing old entry $actual_out");
            unlink $actual_out;
        }
        local *STDOUT = Genome::Sys->open_file_for_writing($actual_out);
        local *STDERR = *STDOUT;
        my @argv = ("model", "define", $sub_command, "-h");

        # using private method here because there isn't a public one that
        # doesn't exit but it is best not to fork or subshell or most
        # specifically to call exit because then we can't track test
        # dependencies
        my $exit = Genome::Command->_cmdline_run(@argv);
        is($exit, 0, 'exited zero: `genome  ' . join(' ', @argv) . '`');
    };

    ok(-s $actual_dir, "output data was generated for $sub_command");

    ok(-s $expected_out, "reference data is available for $sub_command under $expected_dir")
        or next;

    next unless -e $actual_dir;

    my @diff = `diff -r $expected_out $actual_out`;
    #print "diff -r $expected_out $actual_out\n"; <STDIN>;
    is(scalar(@diff), 0, "no differences between actual output and expected output for $sub_command")
        or do {
            for(@diff) { diag($_) };
        }

}

sub module_from_sub_command {
    my $sub_command = shift;

    (my $module = $sub_command) =~ s/-([a-z])/uc($1)/ge;
    $module =~ s/^([a-z])/uc($1)/ge;
    $module = "Genome::Model::$module";

    return $module;
}
