#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More;
require File::Compare;

use_ok( 'Genome::Model::Tools::PooledBac::AddLinkingContigs' ) or die;

my $version = 1;
my $test_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-PooledBac/AddLinkingContigs_v'.$version;
ok( -d $test_dir, 'Test suite dir exist' ) or die;

my $tmp_dir = Genome::Sys->create_temp_directory();
ok( -d $tmp_dir, 'Created temp test dir' ) or die;

my @files_to_copy = qw/
CONTIG_MAP
bac_region_db.blast
bac_region_db.xnd
bac_region_db.xns
bac_region_db.xnt
params.110906.1135
pooled_contigs.fasta
ref_seq.fasta
ref_seq.fasta.qual
/;

for my $file ( @files_to_copy ) {
    ok( -s $test_dir."/$file", "Test $file exists" );
    ok( File::Copy::copy( $test_dir."/$file", $tmp_dir ), "Copied $file to tmp dir" );
    ok( -s $tmp_dir."/$file", "$file exists in tmp dir" );
}

my $tool = Genome::Model::Tools::PooledBac::AddLinkingContigs->create( project_dir => $tmp_dir );
ok( $tool, 'Created tool' );
ok( $tool->execute, 'Successfully executed tool' );

#compare output files
ok( File::Compare::compare($tmp_dir.'/CONTIG_MAP',$test_dir.'/CONTIG_MAP') == 0, 'CONTIG_MAP files match' );
my $report_file = '/reports/contigs_that_link_to_matching_contigs';
ok( File::Compare::compare($tmp_dir."/$report_file",$test_dir."/$report_file")==0, 'Report files match' );

#<STDIN>;

done_testing();
