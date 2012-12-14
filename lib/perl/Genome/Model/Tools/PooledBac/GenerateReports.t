#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More;
require File::Compare;

use_ok( 'Genome::Model::Tools::PooledBac::GenerateReports' ) or die;

my $version = 1;
my $test_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-PooledBac/GenerateReports_v'.$version;
ok( -d $test_dir, 'Test suite dir exists' );

my $tmp_dir = Genome::Sys->create_temp_directory();
ok( -d $tmp_dir, 'Tmp test dir created' );

my @files_to_link = qw/
CONTIG_MAP
H_GD-332K02
bac_region_db.blast
bac_region_db.xnd
bac_region_db.xns
bac_region_db.xnt
params.110906.1135
pooled_contigs.fasta
ref_seq.fasta
ref_seq.fasta.qual
/;

for my $file ( @files_to_link ) {
    ok( -s $test_dir."/$file", "Test dir $file exists" );
    symlink( $test_dir."/$file",$tmp_dir."/$file" );
    ok( -l $tmp_dir."/$file", "Linked $file in temp dir" );
}

my $tool = Genome::Model::Tools::PooledBac::GenerateReports->create( project_dir => $tmp_dir );
ok( $tool, 'Created tool' );
ok( $tool->execute, 'Executed tool' );

#compare output files
my @files_to_compare = qw/
ambiguous_matching_contigs
blast_report
complete_contig_list
complete_contig_list_with_orphan_contigs
contigs_with_multiple_hits
matching_contigs
orphan_contigs
/;

for my $file ( @files_to_compare ) {
    ok( -e $test_dir."/reports/$file", "Test $file exists" );
    ok( -e $tmp_dir."/reports/$file", "New $file created" );
    ok( File::Compare::compare($test_dir."/reports/$file",$tmp_dir."/reports/$file") ==0, "$file files match" );
}

#<STDIN>;

done_testing();
