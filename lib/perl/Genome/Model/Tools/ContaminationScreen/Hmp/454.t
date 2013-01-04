#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More tests => 3;

BEGIN {use_ok('Genome::Model::Tools::ContaminationScreen::Hmp::454');}

my %params;
$params{input_file} = '/gsc/var/tmp/fasta/Hmp/short_454.fna';
$params{database} = '/gscmnt/sata156/research/mmitreva/databases/human_build36/Homo_sapiens.NCBI36.45.dna.aml.plus5.8s15s28s.fna';
$params{output_file} = Genome::Sys->create_temp_file_path('454.screened');
$params{parsed_file} = Genome::Sys->create_temp_file_path('454.parsed');
$params{filter_list} = Genome::Sys->create_temp_file_path('filtered.txt');
$params{percent} = 90;
$params{query_length} = 50;

my $hcs_454 = Genome::Model::Tools::ContaminationScreen::Hmp::454->create(%params);

isa_ok($hcs_454, 'Genome::Model::Tools::ContaminationScreen::Hmp::454');

ok($hcs_454->execute,"454 executing");
