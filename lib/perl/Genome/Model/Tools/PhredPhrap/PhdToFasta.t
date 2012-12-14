#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use File::Compare;
use File::Temp;
use Test::More;

use_ok('Genome::Model::Tools::PhredPhrap::PhdToFasta') or die;

my $path = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-PhredPhrap'; #directory for sample data and output
my $static = "$path/PhdToFasta/_fasta_output.static"; #static file for comparison to output
my $tmpdir = File::Temp::tempdir(CLEANUP => 1);
my $fasta_file = $tmpdir.'/fasta';

my $phd_to_fasta = Genome::Model::Tools::PhredPhrap::PhdToFasta->create(
    phd_dir => "$path/phd_dir/",
    phd_file => "$path/PhdToFasta/phd.txt",  
    fasta_file => $fasta_file,
    _error_file => "$path/PhdToFasta/error.txt", 
);
isa_ok($phd_to_fasta, "Genome::Model::Tools::PhredPhrap::PhdToFasta");
ok($phd_to_fasta->execute,'execute PhdToFasta');
is(File::Compare::compare($phd_to_fasta->fasta_file, "$path/PhdToFasta/_fasta_output.static"), 0, 'fasta file matches');
is(File::Compare::compare($phd_to_fasta->qual_file, "$path/PhdToFasta/_fasta_output.qual.static"), 0, 'qual file matches');

done_testing();
