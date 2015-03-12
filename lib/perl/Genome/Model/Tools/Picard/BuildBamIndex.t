#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More;
use File::Spec qw();

my $pkg = 'Genome::Model::Tools::Picard::BuildBamIndex';
use_ok($pkg);

my $data_dir = File::Spec->catfile($ENV{GENOME_TEST_INPUTS}, 'Genome-Model-Tools-Picard-BuildBamIndex');

my $bam_file = File::Spec->catfile($data_dir, "coordsort.bam");
my $tmpdir = Genome::Sys->create_temp_directory;

my $input_file = File::Spec->catfile($tmpdir, "input.bam");
symlink $bam_file, $input_file;

subtest "default output filename" => sub {
    my $cmd = $pkg->create(input_file => $input_file);
    ok($cmd->execute, "Executed command");
    ok(-s File::Spec->catfile($tmpdir, "input.bai"), "Index created in expected location");
};

subtest "specific output filename" => sub {
    my $output_file = File::Spec->catfile($tmpdir, "foo.bai");
    my $cmd = $pkg->create(input_file => $input_file, output_file => $output_file);
    ok($cmd->execute, "Executed command");
    ok(-s $output_file, "Index created in expected location");
};

done_testing();
