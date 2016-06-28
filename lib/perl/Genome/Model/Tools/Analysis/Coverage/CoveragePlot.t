#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";
use Genome::Utility::Test qw(compare_ok);
use Test::More tests => 18;

eval
{

    # creates a class instance
    my $class = 'Genome::Model::Tools::Analysis::Coverage::CoveragePlot';
    use_ok($class);


    # checks test data files
    my $data_dir = Genome::Utility::Test->data_dir($class);
    ok(-d $data_dir, "data_dir exists: $data_dir") or die;


    # checks inputs
    my $roi_file_path = "$data_dir/test-roi425.bed";
    ok(-s $roi_file_path, 'input exists: roi_file_path') or die;

    my $bam_file1 = "$data_dir/test-bam-CKR.bam";
    ok(-s $bam_file1, 'input exists: bam_files') or die;

    my $bam_file2 = "$data_dir/test-bam-IKR.bam";
    ok(-s $bam_file2, 'input exists: bam_files') or die;

    my $bam_file3 = "$data_dir/test-bam-VKR.bam";
    ok(-s $bam_file3, 'input exists: bam_files') or die;


    #create temp directory for merging
    my $tempdir = Genome::Sys->create_temp_directory();
    ok(-d $tempdir, 'temporary directory created') or die;

    # creates and executes a command
    my $cmd = $class->create(
        roi_file_path => $roi_file_path,
        output_directory => $tempdir,

        bam_files => [$bam_file1, $bam_file2, $bam_file3],
        labels => ["sample1", "sample2", "sample3"],
        min_depth_filters => [1, 20, 40, 60, 80, 100],
        wingspan  => 0,
        band_width => 10,
    );

    ok($cmd, 'created command') or die;
    ok($cmd->execute, 'executed command') or die;


    # check outputs
    ok(-s "$tempdir/table.tsv", 'table.tsv has size');
    compare_ok("$tempdir/table.tsv", "$data_dir/table.tsv", 'table matches');

    ok(-s "$tempdir/boxplot.tab", 'boxplot.tab has size');
    compare_ok("$tempdir/boxplot.tab", "$data_dir/boxplot.tab", 'boxplot matches');

    ok(-s "$tempdir/coverage.tab", 'coverage.tab has size');
    compare_ok("$tempdir/coverage.tab", "$data_dir/coverage.tab", 'coverage matches');

    ok(-s "$tempdir/density.tab", 'density.tab has size');
    compare_ok("$tempdir/density.tab", "$data_dir/density.tab", 'density matches');

    # skips this plot PDF file
    ok(-s "$tempdir/coverage-plot.pdf", 'coverage-plot.pdf has size');
    #compare_ok("$tempdir/coverage-plot.pdf", "$data_dir/coverage-plot.pdf", 'coverage-plot matches');
};

done_testing;
