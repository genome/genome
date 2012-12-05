#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";
use Test::More;
require File::Compare;

use_ok( 'Genome::Model::Tools::Velvet::CreateGapFile' ) or die;

my $version = 'v1';

my $data_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Velvet/CmtFile-'.$version;
ok( -d $data_dir, 'Data dir' );

my $test_dir = Genome::Sys->create_temp_directory();
Genome::Sys->create_directory( $test_dir.'/edit_dir' );
ok( -d $test_dir.'/edit_dir', 'Test dir');

# check/copy test files
my $contigs_bases_file = $data_dir.'/edit_dir/contigs.bases';
ok( -s $contigs_bases_file, 'contigs.bases file' );   
ok( File::Copy::copy( $contigs_bases_file, $test_dir.'/edit_dir' ), 'Copied contigs.bases file' );

my $input_fastq_file = $data_dir.'/input-collated.fastq';
ok( -s $input_fastq_file, 'input fastq file' );
ok( File::Copy::copy( $input_fastq_file, $test_dir ), 'Copied input-collated.fastq file' );

# run tool
my $create = Genome::Model::Tools::Velvet::CmtFile->create(
    version => '1.1.06',
    sequencing_technologies => ['illumina'],
    assembly_directory => $test_dir,
);
ok( $create, 'Created tool' );
ok( $create->execute, 'Executed tool' );

# check output
ok( -s $data_dir.'/edit_dir/contigs.cmt', 'Data dir contigs.cmt file exists' );
ok( -s $test_dir.'/edit_dir/contigs.cmt', 'New contigs.cmt file created' );
ok( File::Compare::compare( $data_dir.'/edit_dir/contigs.cmt', $test_dir.'/edit_dir/contigs.cmt' ) == 0, 'Files match' );

#<STDIN>;

done_testing();

exit;
