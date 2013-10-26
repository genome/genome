#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";
use File::Temp 'tempdir';
use Test::More tests => 7;

BEGIN {
        use_ok('Genome::Model::Tools::Blat::Subject');
}

my $query_file = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Blat-Subject/test.fa';
ok(-s $query_file,'query file has size');

my $subject_file = Genome::Config::reference_sequence_directory() . '/refseq-for-test/11.fa';
ok(-s $subject_file,'subject file has size');

my $output_directory = tempdir(
    'Genome-Model-Tools-Blat-XXXXXX',
    CLEANUP => 1,
    UNLINK => 1,
    TMPDIR => 1,
);

print $output_directory . "\n";
my $blat = Genome::Model::Tools::Blat::Subject->create(
    query_file => $query_file,
    subject_file => $subject_file,
    output_directory => $output_directory,
);

isa_ok($blat,'Genome::Model::Tools::Blat::Subject');
ok($blat->execute,'execute command '. $blat->command_name);
ok($blat->alignment_file =~ /\/test_11\.psl$/,'expected alignment file');
ok(-s $blat->alignment_file,'alignment file has size');
