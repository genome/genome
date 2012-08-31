#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";
use Test::More tests => 3;
use File::Compare;
   
BEGIN
{
    use_ok ('Genome::Model::Tools::Fastq::RemoveN');
}

my $path = $ENV{GENOME_TEST_INPUTS} . "/Genome-Model-Tools-Fastq-RemoveN/";
my $fastq_file      = "$path/in.fastq";
my $n_static        = "$path/expected-out.fastq";
my $n_removed_file  = Genome::Sys->create_temp_file_path('actual-out.fastq');

my $cmd = Genome::Model::Tools::Fastq::RemoveN->create(
    fastq_file     =>  $fastq_file,
    n_removed_file =>  $n_removed_file,
    n_removal_threshold => 1,
);
ok($cmd->execute, "executed without errors");

is(compare($n_removed_file,$n_static),0, "file content from converting $fastq_file to $n_removed_file matches $n_static")
    or do {
    diag("diff $n_static $n_removed_file | head:\n" . `diff $n_static $n_removed_file | head`);
    };


