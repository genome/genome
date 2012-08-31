#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Test::More tests => 3;

use_ok('Genome::Model::Tools::Kmer::Suffixerator');
my $tmp_dir = File::Temp::tempdir('Kmer-Suffixerator-'. Genome::Sys->username .'-XXXX', DIR => "$ENV{GENOME_TEST_TEMP}", CLEANUP => 1);
my $data_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Fasta-ToTwoBit';
my @fasta_files;
for ( 11 .. 13) {
    push @fasta_files, $data_dir .'/'. $_ .'.fasta';
}
my $suffixerator = Genome::Model::Tools::Kmer::Suffixerator->create(
   fasta_files => \@fasta_files,
   index_name=> $tmp_dir .'/test',
   log_file => $tmp_dir .'/test.log',
);
isa_ok($suffixerator,'Genome::Model::Tools::Kmer::Suffixerator');
ok($suffixerator->execute,'execute command '. $suffixerator->command_name);

exit;
