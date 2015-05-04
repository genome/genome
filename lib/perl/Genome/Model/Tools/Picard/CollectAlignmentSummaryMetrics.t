#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More;
use Genome::Utility::Test qw(compare_ok);
use File::Spec qw();

my $pkg = 'Genome::Model::Tools::Picard::CollectAlignmentSummaryMetrics';

use_ok($pkg);

my $data_dir = File::Spec->catfile(Genome::Config::get('test_inputs'), 'Genome-Model-Tools-Picard-CollectAlignmentSummaryMetrics');
my $bam_file = File::Spec->catfile($data_dir, 'aligned.bam');
my $ref_file = File::Spec->catfile($data_dir, 'small.fa');
my $expected_file = File::Spec->catfile($data_dir, 'expected.txt');
my $output_file = Genome::Sys->create_temp_file_path;

my $obj = $pkg->create(
    input_file => $bam_file,
    output_file => $output_file,
    refseq_file => $ref_file,
    );

ok($obj, 'created command');
ok($obj->execute, 'executed command');

compare_ok($expected_file, $output_file, diag => 1, filters => ['^#.*']);
done_testing();
