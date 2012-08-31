#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::More;
use above "Genome";
use File::Basename;
require File::Compare;

use_ok( 'Genome::Model::Tools::Newbler::CreatePlacedReadsFiles' );

my $version = 'v3';
my $test_suite = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Newbler/CreatePlacedReadsFiles-'.$version;
ok( -d $test_suite, "Test suite dir exists" );

my $temp_dir = Genome::Sys->create_temp_directory();
ok( -d $temp_dir, "Temp dir created" );

for my $file ( qw/ 454AllContigs.fna 454ReadStatus.txt 454Scaffolds.txt / ) {
    ok( -s $test_suite."/$file", "Test suite $file file exists" );
    symlink( $test_suite."/$file", $temp_dir."/$file" );
    ok( -l $temp_dir."/$file", "$file file linked in temp dir" );
}

my $tool = Genome::Model::Tools::Newbler::CreatePlacedReadsFiles->create(
    assembly_directory => $temp_dir,
    min_contig_length => 200,
);
ok( $tool, "Created tool" );
ok( $tool->execute, "Executed tool" );

for my $file ( '/consed/edit_dir/readinfo.txt', '/consed/edit_dir/reads.placed' ) {
    ok( -s $test_suite."/$file", "Test suite $file file exists" );
    ok( -s $temp_dir."/$file", "Test created $file file" );
    ok( File::Compare::compare($test_suite."/$file",$temp_dir."/$file")==0, "$file files match" );
}

#<STDIN>;

done_testing();

exit;
