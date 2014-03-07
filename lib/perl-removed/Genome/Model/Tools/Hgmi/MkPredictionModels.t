#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";
use File::Remove qw/ remove /;
use Test::More;

my $d2014_03_04 = 1393975522;
my $d2014_03_11 = 60 * 60 * 24 * 7 + $d2014_03_04;
if (time > $d2014_03_11) {
    plan tests => 8;
} else {
    plan skip_all => 'disabled for phase 3 disk maintenance';
}

use File::Temp 'tempdir';
use Cwd;

use_ok('Genome::Model::Tools::Hgmi::MkPredictionModels');

my $fasta =  $ENV{GENOME_TEST_INPUTS} . "/Genome-Model-Tools-Hgmi/BIFCATDFT.v1.contigs.newname.fasta";
ok(-e $fasta, "fasta file exists at $fasta");

my $test_dir = tempdir(
    'Genome-Model-Tools-Hgmi-MkPredictionModels-XXXXXX',
    TMPDIR => 1,
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
