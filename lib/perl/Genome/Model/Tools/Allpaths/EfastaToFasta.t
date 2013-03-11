#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";

require File::Compare;
use Test::More;

use_ok('Genome::Model::Tools::Allpaths::EfastaToFasta') or die;

my $test_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Allpaths/EfastaToFasta/v1';
my $example_clean_fasta = $test_dir.'/clean.final.assembly.fasta';      

my $temp_dir = Genome::Sys->create_temp_directory();
my $new_clean_fasta = $temp_dir.'/clean.final.assembly.fasta';

my $convert = Genome::Model::Tools::Allpaths::EfastaToFasta->create(
   assembly_directory => $test_dir,
   clean_fasta_out    => $new_clean_fasta,
);

ok( $convert, "create" ) or die;
$convert->dump_status_messages(1);
ok( $convert->execute, "execute" ) or die;

is(File::Compare::compare( $example_clean_fasta, $new_clean_fasta), 0, "files match" );

#print $new_clean_fasta."\n"; <STDIN>;

done_testing();
