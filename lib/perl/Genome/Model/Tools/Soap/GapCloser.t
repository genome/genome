#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More;
require File::Compare;

use_ok( 'Genome::Model::Tools::Soap::GapCloser' ) or die;

#test suite dir
my $version = 2;
my $test_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Soap/GapCloser_v'.$version;
ok( -d $test_dir, 'test dir exists' );

#input fastq files
my @data_files = qw/ 
SRR038746.pair1.fastq
SRR038746.pair2.fastq
SRR038746.single.fastq
SRR042027.pair1.fastq
SRR042027.pair2.fastq
SRR042027.single.fastq
/;
for my $file ( @data_files ) {
    ok( -s $test_dir."/input/$file", "Input $file exists" );
}

#assembly files
my @assembly_files = qw/
config_file
TEST.scafSeq
gapfill
gapfill.fill
/;
for my $file ( @assembly_files ) {
    ok ( -s $test_dir."/$file", "Assembly $file exists" );
}

#make temp directory
my $temp_dir = Genome::Sys->create_temp_directory();
ok( -d $temp_dir, 'temp dir created' );
Genome::Sys->create_directory( $temp_dir.'/edit_dir' );
ok( -d $temp_dir.'/edit_dir', 'Made edit_dir in temp dir' );

#copy files to temp directory
for my $file ( qw/ config_file TEST.scafSeq / ) {
    ok( File::Copy::copy($test_dir."/$file",$temp_dir), "Copied $file to temp dir" );
}

my $fail_2 = Genome::Model::Tools::Soap::GapCloser->create(
    assembly_directory => $temp_dir,
    version => '1.10',
    p => 32,
);
ok( ! $fail_2->execute, 'Failed execute with p > 31' ); 

my $fail_3 = Genome::Model::Tools::Soap::GapCloser->create(
    assembly_directory => $temp_dir,
    p => 25,
);
ok( ! $fail_3->execute, 'Failed execute without version specified' );

#PASS test create/execute tool
my $tool = Genome::Model::Tools::Soap::GapCloser->create(
    assembly_directory => $temp_dir,
    version => '1.10',
    p => 25,
);
$tool->dump_status_messages(1);
ok( $tool, "Created tool" );
ok( $tool->execute, "Executed tool" );

#check/compare output
for my $file ( qw/ gapfill gapfill.fill / ) {
    ok( -s $temp_dir."/edit_dir/$file", "Created $file" );
    ok( File::Compare::compare($temp_dir."/edit_dir/$file",$test_dir."/$file")==0, "$file files match" );
}

#<STDIN>;

done_testing();
