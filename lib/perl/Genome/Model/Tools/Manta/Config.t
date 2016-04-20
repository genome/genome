#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More tests => 3;

my $pkg = 'Genome::Model::Tools::Manta::Config';
use_ok ($pkg);

my $version = '0.29.6';
my $expected_cmdline = '/gscmnt/gc13001/info/model_data/jwalker_scratch/src/manta-'. $version .'.centos5_x86_64/bin/configManta.py --bam "/tmp/bam1" --bam "/tmp/bam2" --bam "/tmp/bam3" --referenceFasta "/tmp/ref.fa"';

my @bams = ('/tmp/bam1','/tmp/bam2','/tmp/bam3');
my $fasta = '/tmp/ref.fa';
my @expected_input_files = (@bams,$fasta);
my $config_manta = Genome::Model::Tools::Manta::Config->create(
   version => $version,
   bams => \@bams,
   reference_fasta => $fasta,
);

my $cmdline_string = $config_manta->build_cmdline_string();
is($cmdline_string,$expected_cmdline,'expected command-line for version');

$config_manta->_resolve_input_files;
my @input_files = $config_manta->input_files;                       
is_deeply(\@input_files,\@expected_input_files,'expected input files');

# TODO: Test the execute subroutine

exit;