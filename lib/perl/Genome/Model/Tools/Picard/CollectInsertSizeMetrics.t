#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More;
use Genome::Utility::Test qw(compare_ok);
use File::Spec qw();

my $pkg = 'Genome::Model::Tools::Picard::CollectInsertSizeMetrics';

use_ok($pkg);

my $data_dir = File::Spec->catfile($ENV{GENOME_TEST_INPUTS}, 'Genome-Model-Tools-Picard-CollectInsertSizeMetrics');
my $bam_file = File::Spec->catfile($data_dir, 'aligned.bam');
my $expected_file = File::Spec->catfile($data_dir, 'expected.txt');
my $output_file = Genome::Sys->create_temp_file_path;
my $plot_file = Genome::Sys->create_temp_file_path;

my $obj = $pkg->create(
    input_file => $bam_file,
    output_file => $output_file,
    histogram_file => $plot_file,
    # ^^^^ don't diff this. different R versions/backends may differ
    # in the exact details of how they write pdf files.
    # (it's smart to diff graphics files. no, wait, the other thing)
    );

ok($obj, 'created command');
ok($obj->execute, 'executed command');

compare_ok($expected_file, $output_file, diag => 1, filters => ['^#.*']);
done_testing();
