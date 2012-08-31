#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More tests => 2;
use FindBin qw($Bin);

my $datadir = $Bin . '/454.t.d';

BEGIN {use_ok('Genome::Model::Tools::ContaminationScreen::454');}

my %params;
$params{input_file} = $datadir . '/test.fna';
$params{database} = '/gscmnt/sata156/research/mmitreva/databases/human_build36/HS36.chr_Mt_ribo.fna';
my $hcs_454 = Genome::Model::Tools::ContaminationScreen::454->create(%params);

isa_ok($hcs_454, 'Genome::Model::Tools::ContaminationScreen::454');




