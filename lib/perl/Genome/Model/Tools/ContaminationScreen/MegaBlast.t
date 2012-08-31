#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More tests => 2;


BEGIN {use_ok('Genome::Model::Tools::ContaminationScreen::MegaBlast');}

my %params;
$params{input_file} = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-ContaminationScreen-MegaBlast/test_nt.fna';
$params{output_file} = Genome::Sys->create_temp_file_path('test_output.fna');
$params{database} = '/gscmnt/sata837/assembly/nt_db/genbank_nt_20091004'; 
$params{header} = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-ContaminationScreen-MegaBlast/nt.index.header';

my $hcs_MegaBlast = Genome::Model::Tools::ContaminationScreen::MegaBlast->create(%params);

isa_ok($hcs_MegaBlast, 'Genome::Model::Tools::ContaminationScreen::MegaBlast');
