#!/usr/bin/env genome-perl

use strict;
use warnings;

use File::Compare;
use Test::More tests => 4;

use above 'Genome';

use_ok('Genome::Model::Tools::Kmer::DuplicationRate');

my $tmp_dir = Genome::Sys->create_temp_directory('Genome-Model-Tools-Kmer-PlotDuplicationRate-'. Genome::Sys->username);
my $data_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Kmer-DuplicationRate';
my $occratio_file = $data_dir .'/s_7_occratio.txt';

my $plot_dup_rate = Genome::Model::Tools::Kmer::PlotDuplicationRate->create(
    occratio_file => $occratio_file,
    plot_file => $tmp_dir .'/s_7_occratio.png',
);
isa_ok($plot_dup_rate,'Genome::Model::Tools::Kmer::PlotDuplicationRate');
ok($plot_dup_rate->execute,'execute command '. $plot_dup_rate->command_name);
ok(-s $plot_dup_rate->plot_file,'Found plot file '. $plot_dup_rate->plot_file);
