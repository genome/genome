#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";
use Genome::Utility::Test qw(compare_ok);
use Test::More tests => 11;

eval
{
    
    # creates a class instance
    my $class = 'Genome::Model::Tools::Analysis::Coverage::MergeReadcounts';
    use_ok($class);
    
    
    # checks test data files
    my $data_dir = Genome::Utility::Test->data_dir($class);
    # $data_dir = './3.test_suite_data';     # temporarily switched for debug
    ok(-d $data_dir, "data_dir exists: $data_dir") or die;
    
    
    # checks inputs
    my $bam_file1 = "$data_dir/merge-readcounts-sample1-chr1.bam";
    ok(-s $bam_file1, 'input exists: bam_file1') or die;
    my $bam_file2 = "$data_dir/merge-readcounts-sample2-chr1.bam";
    ok(-s $bam_file2, 'input exists: bam_file2') or die;
    
    my $variant_file1 = "$data_dir/merge-readcounts-snvs.indels.example1-chr1";
    ok(-s $variant_file1, 'input exists: variant_file1') or die;
    my $variant_file2 = "$data_dir/merge-readcounts-snvs.indels.example2-chr1";
    ok(-s $variant_file2, 'input exists: variant_file2') or die;
    
    my $output_file = "$data_dir/output-snvs.indels.annotated-merge";
    ok(-s $output_file, 'input exists: output_file') or die;
    
    
    #create temp directory for merging
    my $tempdir = Genome::Sys->create_temp_directory();
    #ok(-s $tempdir, 'temporary directory created') or die;
    
    # creates and executes a command
    my $cmd = $class->create(
        bam_files => "$bam_file1,$bam_file2",
        variant_files => "$variant_file1,$variant_file2",
        variant_sources => "sample1,sample2",
        genome_build => "mm9",
        output_file => "$tempdir/out_file",
        header_prefixes => "normal,tumor",
    );
    
    ok($cmd, 'created command') or die;
    ok($cmd->execute, 'executed command') or die;
    
    
    # check outputs
    ok(-s "$tempdir/out_file", 'output_file has size');
    compare_ok("$tempdir/out_file", $output_file, 'output matches');
};

done_testing;
