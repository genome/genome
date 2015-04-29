#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More tests => 2;

BEGIN {use_ok('Genome::Model::Tools::ContaminationScreen::Hmp::Illumina');}

my %params;
$params{dir} = Genome::Config::get('test_inputs') . '/Genome-Model-Tools-ContaminationScreen-Hmp-Illumina';#
$params{fastq1} = Genome::Config::get('test_inputs') . '/Genome-Model-Tools-ContaminationScreen-Hmp-Illumina/solexa-fastq-contam0'; 
$params{fastq2} = Genome::Config::get('test_inputs') . '/Genome-Model-Tools-ContaminationScreen-Hmp-Illumina/solexa-fastq-contam1';

my $illumina = Genome::Model::Tools::ContaminationScreen::Hmp::Illumina->create(%params);

isa_ok($illumina, 'Genome::Model::Tools::ContaminationScreen::Hmp::Illumina');

#$illumina->execute();
