#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";
use Test::More tests => 3;
use File::Compare;
use File::Temp;

BEGIN {
        use_ok ('Genome::Model::Tools::Fastq::ToPhdballScf');
}

my $dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Fastq/ToPhdballScf';
my $fastq_file = $dir .'/test.fq';
my $time = 'Mon Jan 10 20:00:00 2009';

my $tmp_dir = File::Temp::tempdir(
    "PhdBallScf_XXXXXX", 
    TMPDIR => 1,
    CLEANUP => 1,
);

my $scf_dir   = $tmp_dir.'/chromat_dir';
my $ball_file = $tmp_dir.'/phd.ball';

my %params = (
    fastq_file => $fastq_file,
    scf_dir    => $scf_dir,
    ball_file  => $ball_file,
    id_range   => '1-33',
    base_fix   => 1,
    time       => $time,
    solexa_fastq => 1,
);

my $to_ballscf = Genome::Model::Tools::Fastq::ToPhdballScf->create(%params);

isa_ok($to_ballscf,'Genome::Model::Tools::Fastq::ToPhdballScf');
ok($to_ballscf->execute,'ToPhdballScf executes ok');
