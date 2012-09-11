#!/gsc/bin/perl

use strict;
use warnings;

use IO::File;
use above 'Genome';

use Test::More;

plan tests=>8;

my $test_dir = $ENV{GENOME_TEST_INPUTS} . "/Genome-Model-Tools-Relationship-BladeRunner/";
my $job1_stdout = $test_dir . "/job1_std_out";
my $job1_stderr = $test_dir . "/job1_std_err";
my $job2_stdout = $test_dir . "/job2_std_out";
my $job2_stderr = $test_dir . "/job2_std_err";

my $job2 = Genome::Model::Tools::Relationship::LsfJob->create(command_line=>"rm $test_dir/touch_testfile", std_out=>"$job2_stdout", std_err=>"$job2_stderr");
my $job = Genome::Model::Tools::Relationship::LsfJob->create(command_line=>"touch /gscuser/charris/bladerunner_test", std_out=>"$job1_stdout", std_err=>"$job1_stderr", parents=>[$job2]);
my @jobs = ($job, $job2);

ok(!eval{ Genome::Model::Tools::Relationship::BladeRunner->execute(jobs=>\@jobs)}, "Blade Runner Failed-- and it should have");
ok(-e $job2_stderr, "First jobs error output exists!");
ok(-e $job2_stdout, "First jobs standard output exists!");
`rm -f $job2_stderr $job2_stdout`;

$job = Genome::Model::Tools::Relationship::LsfJob->create(command_line=>"touch $test_dir/touch_testfile", std_out=>"$job1_stdout", std_err=>"$job1_stderr");
$job2 = Genome::Model::Tools::Relationship::LsfJob->create(command_line=>"rm $test_dir/touch_testfile", std_out=>"$job2_stdout", std_err=>"$job2_stderr", parents=>[$job]);
@jobs = ($job, $job2);
ok(Genome::Model::Tools::Relationship::BladeRunner->execute(jobs=>\@jobs), "Blade Runner Succeeded on a simple 2job workflow");
ok(-e $job2_stderr, "Second jobs error output exists!");
ok(-e $job2_stdout, "Second jobs standard output exists!");
ok(-e $job1_stderr, "First jobs error output exists!");
ok(-e $job1_stdout, "First jobs standard output exists!");
`rm -f $job2_stderr $job2_stdout $job1_stdout $job1_stderr`;



