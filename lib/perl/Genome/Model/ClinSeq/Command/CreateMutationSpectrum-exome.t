#!/usr/bin/env genome-perl

#Written by Malachi Griffith

use strict;
use warnings;
use File::Basename;
use Cwd 'abs_path';

BEGIN {
    $ENV{UR_DBI_NO_COMMIT}               = 1;
    $ENV{UR_USE_DUMMY_AUTOGENERATED_IDS} = 1;
}

use above "Genome";
use Test::More tests => 11;  #One per 'ok', 'is', etc. statement below
use Genome::Model::ClinSeq::Command::CreateMutationSpectrum;
use Data::Dumper;

use_ok('Genome::Model::ClinSeq::Command::CreateMutationSpectrum') or die;

#Define the test where expected results are stored
my $expected_output_dir =
    Genome::Config::get('test_inputs') . "/Genome-Model-ClinSeq-Command-CreateMutationSpectrum/exome/2014-10-10/";
ok(-e $expected_output_dir, "Found test dir: $expected_output_dir") or die;

#Create a temp dir for results
my $temp_dir = Genome::Sys->create_temp_directory();
ok($temp_dir, "created temp directory: $temp_dir");

#Get a wgs somatic variation build
my $clinseq_build_id = "27b94a6da5a44520bf7816ac26650f6e";
my $clinseq_build    = Genome::Model::Build->get($clinseq_build_id);
ok($clinseq_build, "obtained clinseq build from db") or die;

#Get a wgs somatic variation build
my $somvar_build_id = 129399487;
my $somvar_build    = Genome::Model::Build->get($somvar_build_id);
ok($somvar_build, "obtained somatic variation build from db") or die;

#Get a 'final' name for the sample
my $final_name = $somvar_build->model->id;
$final_name = $somvar_build->model->subject->name if ($somvar_build->model->subject->name);
$final_name = $somvar_build->model->subject->individual->common_name
    if ($somvar_build->model->subject->individual->common_name);
ok($final_name, "got final name from build object") or die;

#Create create-mutation-spectrum command and execute
#genome model clin-seq create-mutation-spectrum --outdir=/tmp/create_mutation_spectrum/ --datatype=exome --max-snvs=100 129346170

my $mutation_spectrum_cmd = Genome::Model::ClinSeq::Command::CreateMutationSpectrum->create(
    outdir        => $temp_dir,
    datatype      => "exome",
    clinseq_build => $clinseq_build,
    somvar_build  => $somvar_build,
    test          => 1
);
$mutation_spectrum_cmd->queue_status_messages(1);
my $r1 = $mutation_spectrum_cmd->execute();
is($r1, 1, 'Testing for successful execution.  Expecting 1.  Got: ' . $r1);

#Dump the output to a log file
my @output1  = $mutation_spectrum_cmd->status_messages();
my $log_file = $temp_dir . "/exome/CreateMutationSpectrum.exome.log.txt";
my $log      = IO::File->new(">$log_file");
$log->print(join("\n", @output1));
ok(-e $log_file, "Wrote message file from update-analysis to a log file: $log_file");

#The first time we run this we will need to save our initial result to diff against
#Genome::Sys->shellcmd(cmd => "cp -r -L $temp_dir/* $expected_output_dir");

#Check for non-zero presence of expected PDFs
my $pdf3 =
      $temp_dir
    . "/exome/mutation_spectrum_sequence_context/"
    . "$final_name"
    . ".mutation-spectrum-sequence-context.pdf";
ok(-s $pdf3, "Found non-zero PDF file mutation-spectrum-sequence-context.pdf");

my $pdf4 = $temp_dir . "/exome/summarize_mutation_spectrum/" . "$final_name" . "_summarize-mutation-spectrum.pdf";
ok(-s $pdf4, "Found non-zero PDF file summarize-mutation-spectrum.pdf");

#Perform a diff between the stored results and those generated by this test
my @diff =
    `diff -r -x '*.log.txt' -x '*.pdf' -x '*.stderr' -x '*.stdout' -x '*.input.tsv' $expected_output_dir $temp_dir`;
ok(@diff == 0, "Found only expected number of differences between expected results and test results")
    or do {
    diag("expected: $expected_output_dir\nactual: $temp_dir\n");
    diag("differences are:");
    diag(@diff);
    my $diff_line_count = scalar(@diff);
    print "\n\nFound $diff_line_count differing lines\n\n";
    Genome::Sys->shellcmd(cmd => "rm -fr /tmp/last-create-mutation-spectrum-exome-result/");
    Genome::Sys->shellcmd(cmd => "mv $temp_dir /tmp/last-create-mutation-spectrum-exome-result");
    };

