#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More; 

use_ok('Genome::Model::Tools::PooledBac::RunBlast') or die;

my $version = 2;
my $test_dir = Genome::Config::get('test_inputs') . '/Genome-Model-Tools-PooledBac/RunBlast_v'.$version;
ok( -d $test_dir, 'Test suite dir exists' ) or die;

my $run_dir = Genome::Sys->create_temp_directory();
ok( -d $run_dir, 'Temp test dir created' );

my %params = (
    ace_file_name => 'ace.dft',
    pooled_bac_dir => $test_dir.'/pooled_bac',
    ref_seq_file => $test_dir.'/ref_seq/ref_small.txt',
    project_dir => $run_dir,
);

my $tool = Genome::Model::Tools::PooledBac::RunBlast->create( %params );
ok( $tool, 'Created tool' );
ok( $tool->execute, 'Executed tool' );

#check blast out files exist .. contents could vary
#so just check that they have been made
my @blast_files = qw/
bac_region_db.blast
bac_region_db.xnd
bac_region_db.xns
bac_region_db.xnt
pooled_contigs.fasta
ref_seq.fasta
ref_seq.fasta.qual
/;

for my $file( @blast_files ) {
    ok( -s $run_dir."/$file", "$file created" );
}

#<STDIN>;

done_testing();
