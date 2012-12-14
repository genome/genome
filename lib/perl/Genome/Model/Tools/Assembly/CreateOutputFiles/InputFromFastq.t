#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";
use Test::More;

use_ok( 'Genome::Model::Tools::Assembly::CreateOutputFiles::InputFromFastq' ) or die;

#check data files
my $data_dir = $ENV{GENOME_TEST_INPUTS} . "/Genome-Model-Tools-Assembly-CreateOutputFiles";
ok(-d $data_dir, "Found data directory: $data_dir") or die;
my $test_fastq_file = $data_dir.'/test.fastq';
ok(-s $test_fastq_file, "Found test contigs.fa file") or die;

#make test dir
my $temp_dir = Genome::Sys->create_temp_directory();
Genome::Sys->create_directory( $temp_dir.'/edit_dir' );

#copy input file
ok(File::Copy::copy($test_fastq_file, $temp_dir),"Copied input contigs file");
my $temp_fastq = $temp_dir.'/test.fastq';
ok(-s $temp_fastq, "New temp fastq file exists");

#run
my $ec = system("chdir $temp_dir; gmt assembly create-output-files input-from-fastq --directory $temp_dir --fastq-file $temp_fastq");
ok($ec == 0, "Command ran successfully");

#fasta files
ok( -s $data_dir.'/edit_dir/test.fasta.gz', "Test fasta.gz exists" );
ok( -s $temp_dir.'/edit_dir/test.fasta.gz', "New test.fasta.gz file created" );
my @diff = `zdiff $data_dir/edit_dir/test.fasta.gz $temp_dir/edit_dir/test.fasta.gz`;
is( scalar @diff, 0, "Fasta files match");

#qual files
ok( -s $data_dir.'/edit_dir/test.fasta.qual.gz', "Test fasta.qual.gz file exists" );
ok( -s $temp_dir.'/edit_dir/test.fasta.qual.gz', "New test.fasta.qual.gz file created" );

done_testing();
