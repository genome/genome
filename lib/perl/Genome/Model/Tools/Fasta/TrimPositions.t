#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";
use Test::More;
require File::Compare;

use_ok('Genome::Model::Tools::Fasta::TrimPositions');

my $data_dir = $ENV{GENOME_TEST_INPUTS}.'/Genome-Model-Tools-Fasta/TrimPositions';
ok( -d $data_dir, "data dir" );

my $test_dir = Genome::Sys->create_temp_directory();
ok( -d $test_dir, "made test dir" );

my @test_files = qw/
trim_positions
contigs.bases
contigs.quals
contigs.bases.clipped
contigs.quals.clipped
/;

for my $file ( @test_files ) {
    ok( -s "$data_dir/$file", "test $file exists" );
}

# copy test files
for my $file ( qw/ contigs.bases contigs.quals trim_positions / ) {
    ok( File::Copy::copy( "$data_dir/$file", $test_dir ), "copied $file" );
}

# run tool
my %cmd_params = (
   fasta_file        => $test_dir.'/contigs.bases',
   qual_file         => $test_dir.'/contigs.quals',
   trim_file         => $test_dir.'/trim_positions',
   min_contig_length => 20,
);

my $tool = Genome::Model::Tools::Fasta::TrimPositions->create( %cmd_params );
ok( $tool, "create tool" );

ok( $tool->execute, "execute tool" );

# check output
for my $file ( qw/ contigs.bases.clipped contigs.quals.clipped / ) {
    ok( -s "$test_dir/$file", "made $file" );
    ok( File::Compare::compare("$data_dir/$file", "$test_dir/$file") == 0, "$file matches" );
}

#print $test_dir."\n"; <STDIN>;

done_testing();
