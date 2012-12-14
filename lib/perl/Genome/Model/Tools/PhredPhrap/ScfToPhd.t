#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use File::Compare;
use File::Temp;
use Test::More;

use_ok('Genome::Model::Tools::PhredPhrap::ScfToPhd') or die;

my $path = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-PhredPhrap'; #directory for sample data and output

my $tmpdir = File::Temp::tempdir(CLEANUP => 1);
my $phd_dir = $tmpdir.'/phd_dir';
mkdir $phd_dir;
ok(-d $phd_dir, 'created phd dir') or die;
my $phd_file = $tmpdir.'/scf.phd';

my $scf_to_phd = Genome::Model::Tools::PhredPhrap::ScfToPhd->create(
    phd_dir => $phd_dir,
    chromat_dir => "$path/chromat_dir/",
    scf_file => "$path/ScfToPhd/scf.txt",
    phd_file => $phd_file, 
);
isa_ok($scf_to_phd,'Genome::Model::Tools::PhredPhrap::ScfToPhd');
ok($scf_to_phd->execute,'execute ScfToPhd');
is(File::Compare::compare($scf_to_phd->phd_file, "$path/ScfToPhd/phd.txt"), 0, 'phds file is the same');

done_testing();
