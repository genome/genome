#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Genome::Utility::Test 'compare_ok';
use Test::More;
use File::Spec qw();

my $pkg = 'Genome::Model::Tools::Picard::ConcatenateMetricsFiles';
use_ok($pkg);

my $data_dir = File::Spec->canonpath( sprintf "%s.d", __FILE__ );

my @metrics = qw/
   CollectAlignmentSummaryMetrics
   CollectInsertSizeMetrics
   CalculateHsMetrics
/;

my @prefix = qw/A B C D/;

for my $metrics (@metrics) {
    my @metrics_files;
    for my $prefix (@prefix) {
        my $metrics_file = File::Spec->catfile($data_dir, $prefix .'.'. $metrics .'.txt');
        push @metrics_files, $metrics_file;
    }    

    my $expected_file = File::Spec->catfile($data_dir, 'expected.'. $metrics .'.txt');
    my $output_file = Genome::Sys->create_temp_file_path;
    
    my $cmd = $pkg->create(
       metrics_files => \@metrics_files,
       output_file => $output_file,
    );

    ok($cmd, "created command");
    ok($cmd->execute, "executed command");

    compare_ok($expected_file, $output_file, replace => [[$data_dir, 'DIRECTORY']]);
}

done_testing();
