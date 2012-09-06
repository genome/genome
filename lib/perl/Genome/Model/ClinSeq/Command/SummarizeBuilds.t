#!/usr/bin/env genome-perl

#Written by Malachi Griffith and Scott Smith

use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
    $ENV{UR_USE_DUMMY_AUTOGENERATED_IDS} = 1;
    $ENV{NO_LSF} = 1;
};

use above "Genome";
use Test::More tests=>6; #One per 'ok' statement below
use Genome::Model::ClinSeq::Command::SummarizeBuilds;
use Data::Dumper;

#Create a temp dir for results
my $temp_dir = Genome::Sys->create_temp_directory();
ok($temp_dir, "created temp directory: $temp_dir");

#Pull an existing clinseq build for test purposes
my $b = Genome::Model::Build->get(119971814);
ok($b, "obtained a clinseq build from the database");

#Execute the tool code
#To make the test run faster, use the 'skip_lims_report' option

my $cmd = Genome::Model::ClinSeq::Command::SummarizeBuilds->create(builds=>[$b], outdir=>$temp_dir, skip_lims_reports=>1);
$cmd->queue_status_messages(1);
my $r = $cmd->execute();
ok($r, "execute summarize-builds successfully");

my @output = $cmd->status_messages();
my $log = IO::File->new(">$temp_dir/status_messages.txt");
ok($log, "created a file to hold the status messages") or die;
$log->print(join("\n", @output));

my $expected_data_directory = $ENV{"GENOME_TEST_INPUTS"} . '/Genome-Model-ClinSeq-Command-SummarizeBuilds/2012-09-06';
ok(-d $expected_data_directory, "found expected data directory: $expected_data_directory");

my @diff = `diff -r $expected_data_directory $temp_dir`;
is(@diff, 5, "no differences from expected results and actual other than the 5 lines corresponding to the random temp dir statement")
  or do { 
      diag("differences are:");
      diag(@diff);
      if (-e "/tmp/last-summarize-builds-test-result"){
        Genome::Sys->shellcmd(cmd => "rm -fr /tmp/last-summarize-builds-test-result");
      }
      Genome::Sys->shellcmd(cmd => "mv $temp_dir /tmp/last-summarize-builds-test-result");
      print "\n\nPlaced new differing results in /tmp/last-summarize-builds-test-result\n";
  };

#The summarize-builds tool changes the current working directory because the LIMS illumina_info tool writes output there and has no option to control this behavior
#We need to change to somewhere else so that the testing harness can clean up the temp file created
chdir("/");

