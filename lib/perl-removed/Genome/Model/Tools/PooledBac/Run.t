#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More;
require File::Compare;

if (Genome::Config->arch_os ne 'x86_64') {
    plan skip_all => 'requires 64-bit machine';
}

use_ok( 'Genome::Model::Tools::PooledBac::Run' ) or die;

my $version = 1;
my $test_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-PooledBac/Run_v'.$version;
ok( -d $test_dir, 'Test suite dir exists' ) or die;

my $run_dir = Genome::Sys->create_temp_directory();
ok( -d $run_dir, 'Test run dir created' );

my %params = (
    ace_file_name => 'ace.dft',
    pooled_bac_dir => $test_dir.'/pooled_bac',
    ref_seq_file => $test_dir.'/ref_seq/ref_small.txt',
    project_dir => $run_dir,
);

my $tool = Genome::Model::Tools::PooledBac::Run->create( %params );
ok( $tool, 'Created tool' ) or die;
ok( $tool->execute, 'Successfully executed tool' ) or die;

#compare output files: ace
my $project_name = 'H_GD-332K02';
my $project_ace = 'H_GD-332K02.fasta.screen.ace';
my $old_ace = $test_dir."/project/$project_name/edit_dir/$project_ace";
ok( -s $old_ace, 'Old ace exists' );
my $new_ace = $run_dir."/$project_name/edit_dir/$project_ace";
ok( -s $new_ace, 'New ace created' );
ok( File::Compare::compare($old_ace, $new_ace)==0, 'Ace files match' );

#compare output files: report
my @report_files = qw/
orphan_contigs
contigs_with_multiple_hits
contigs_that_link_to_matching_contigs
ambiguous_matching_contigs
contigs_only_consensus
assembly_size_report
complete_contig_list_with_orphan_contigs
complete_contig_list
contig_size_report
matching_contigs
matching_contigs
/;

for my $file ( @report_files ) {
    ok( -e $test_dir."/project/reports/$file", "Test dir $file exists" );
    ok( -e $run_dir."/reports/$file", "New $file created" );
    ok( File::Compare::compare($test_dir."/project/reports/$file",$run_dir."/reports/$file")==0, "$file files match" );
}

#<STDIN>;

done_testing();
