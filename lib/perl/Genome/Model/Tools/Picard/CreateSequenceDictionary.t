#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Genome::Utility::Test 'compare_ok';
use Test::More;
use Test::Exception;
use File::Spec qw();

my $pkg = 'Genome::Model::Tools::Picard::CreateSequenceDictionary';
use_ok($pkg);

my $data_dir = sprintf "%s.d", __FILE__;

my $ref_file = File::Spec->catfile($data_dir, "ref.fa");
my $expected_file = File::Spec->catfile($data_dir, "expected.sam");
my $output_file = Genome::Sys->create_temp_file_path;

my $cmd = $pkg->create(
    reference_fasta => $ref_file,
    output_file => $output_file
    );

ok($cmd, "created command");
ok($cmd->execute, "executed command");

compare_ok($expected_file, $output_file, filters => ['file:\S*']);

done_testing();
