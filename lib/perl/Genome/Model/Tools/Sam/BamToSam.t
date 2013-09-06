#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::More;
use File::Compare;
use File::Copy;
use File::Temp;

use above 'Genome';

if (`uname -a` =~ /x86_64/){
    plan tests => 3;
} else{
    plan skip_all => 'Must run on a 64 bit machine';
}

use_ok('Genome::Model::Tools::Sam::BamToSam');

my $data_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Sam-BamToSam';

my $tmp_dir  = File::Temp::tempdir(
    "BamToSam_XXXXXX", 
    TMPDIR => 1,
    CLEANUP => 1,
);

my $bam_file = $tmp_dir .'/test.bam';
copy $data_dir.'/test.bam', $bam_file;

my $result = Genome::Model::Tools::Sam::BamToSam->execute(bam_file => $bam_file);

ok(-e $tmp_dir."/test.sam","Found output file, properly named.");
ok($result,"Tool exited properly.");

