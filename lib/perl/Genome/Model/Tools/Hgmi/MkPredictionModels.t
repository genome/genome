#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";
use File::Remove qw/ remove /;
use Test::More tests => 8;
use File::Temp 'tempdir';
use Cwd;

BEGIN {
        use_ok('Genome::Model::Tools::Hgmi::MkPredictionModels');
}

my $fasta =  $ENV{GENOME_TEST_INPUTS} . "/Genome-Model-Tools-Hgmi/BIFCATDFT.v1.contigs.newname.fasta";
ok(-e $fasta, "fasta file exists at $fasta");

my $test_dir = tempdir(
    'Genome-Model-Tools-Hgmi-MkPredictionModels-XXXXXX',
    DIR => "$ENV{GENOME_TEST_TEMP}/",
    UNLINK => 1,
    CLEANUP => 1,
);

my $cwd = getcwd();
chdir($test_dir);

my $m = Genome::Model::Tools::Hgmi::MkPredictionModels->create(
    locus_tag => "BIFCATDFT",
    fasta_file => $fasta,
    work_directory => $test_dir
);

isa_ok($m,'Genome::Model::Tools::Hgmi::MkPredictionModels');
ok($m->execute(),'create models executed successfully');

is($m->gc(), 56, "gc content matches expected value");

my $new_fasta = $ENV{GENOME_TEST_INPUTS} . "/Genome-Model-Tools-Hgmi/BIFCATDFT.v1.contigs.newname.57gc.fasta";
ok(-e $new_fasta, "fasta file exists at $new_fasta");

$m = Genome::Model::Tools::Hgmi::MkPredictionModels->create(
    locus_tag => "BIFCATDFT",
    fasta_file => $new_fasta,
    work_directory => $test_dir,
);

ok($m->execute(),'create models executed successfully!');

is($m->gc(), 59, "gc content matches expected value");

chdir($cwd);
