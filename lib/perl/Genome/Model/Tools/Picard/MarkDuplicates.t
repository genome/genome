#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";
use Genome::Model::Tools::Picard::MarkDuplicates;
use Test::More tests => 6;

my $input = $ENV{GENOME_TEST_INPUTS} . '/Genome-Tools-Sam-MarkDuplicates/sample.bam';

# step 1: test 1 

my $tmp_dir = File::Temp->newdir(
    "MarkDuplicates_XXXXX",
    TMPDIR => 1,
    CLEANUP => 1
);


my $output_file = File::Temp->new(SUFFIX => ".bam", DIR => $tmp_dir);
my $metrics_file = File::Temp->new(SUFFIX => ".metrics", DIR => $tmp_dir);

my $cmd_1 = Genome::Model::Tools::Picard::MarkDuplicates->create(
    input_file => $input,
    output_file => $output_file->filename,
    metrics_file => $metrics_file->filename,
    temp_directory => $tmp_dir->dirname,
    remove_duplicates => 1,
    use_version => "1.85",
);


ok($cmd_1, "created command");
ok($cmd_1->execute, "executed");
ok(-s $output_file->filename, "output file is nonzero");
unlink(map({$_->filename} $output_file, $metrics_file));

my $cmd_2 = Genome::Model::Tools::Picard::MarkDuplicates->create(
    input_file => $input,
    output_file => $output_file->filename,
    metrics_file => $metrics_file->filename,
    temp_directory => $tmp_dir->dirname,
    remove_duplicates => 1,
    use_version => "1.36",
);

ok($cmd_2, "created command");
ok($cmd_2->execute, "executed");
ok(-s $output_file->filename, "output file is nonzero");
