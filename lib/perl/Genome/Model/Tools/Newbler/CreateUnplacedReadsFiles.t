#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::More;
use above "Genome";
require File::Compare;

use_ok( 'Genome::Model::Tools::Newbler::CreateUnplacedReadsFiles' ) or die;

my $version = 'v1';
my $example_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Newbler/CreateUnplacedReadsFiles-'.$version;
ok( -d $example_dir, 'Example dir exists' ) or die;

my $test_dir = Genome::Sys->create_temp_directory();
ok( -d $test_dir, 'Created temp test dir' ) or die;
Genome::Sys->create_directory( $test_dir.'/consed/edit_dir' );
ok( -d $test_dir.'/consed/edit_dir', 'Created consed edit_dir' ) or die;

my @files_to_link = qw/
    2854709902.partial-input.fastq
    454AllContigs.fna
    454ReadStatus.txt
    454Scaffolds.txt
/;
for my $file ( @files_to_link ) {
    ok( -s "$example_dir/$file", "Example $file file exists" );
    symlink( "$example_dir/$file", "$test_dir/$file" );
    ok( -l "$test_dir/$file", "Linked $file file to temp test dir" );
}

#creat/execute
my $tool = Genome::Model::Tools::Newbler::CreateUnplacedReadsFiles->create(
    assembly_directory => $test_dir,
    min_contig_length => 20,
);
ok( $tool, 'Created tool' );
ok( $tool->execute, 'Executed tool' );

#compare output files
for my $file ( qw/ reads.unplaced reads.unplaced.fasta / ) {
    ok( -s "$example_dir/consed/edit_dir/$file", "Example $file file exists" );
    ok( -s "$test_dir/consed/edit_dir/$file", "Test created $file file" );
    ok( File::Compare::compare("$example_dir/consed/edit_dir/$file", "$test_dir/consed/edit_dir/$file") == 0, "$file files match" );
}

#<STDIN>;

done_testing();
