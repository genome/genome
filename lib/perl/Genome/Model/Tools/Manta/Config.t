#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More tests => 3;

my $pkg = 'Genome::Model::Tools::Manta::Config';
use_ok ($pkg);

my $version = '0.29.6';

my @bams = ('/tmp/bam1.bam','/tmp/bam 2.bam','/tmp/bam3.bam');
my $fasta = '/tmp/ref.fa';
my @expected_input_files = (@bams,$fasta);

my $tmp_dir = Genome::Sys->create_temp_directory();

my @expected_cmdline_array = (
   '/gscmnt/gc13001/info/model_data/jwalker_scratch/src/manta-'. $version .'.centos5_x86_64/bin/configManta.py',
   ('--bam', '/tmp/bam1.bam', '--bam', '/tmp/bam 2.bam', '--bam', '/tmp/bam3.bam', '--referenceFasta', '/tmp/ref.fa', '--runDir', $tmp_dir),
);

my $config_manta = Genome::Model::Tools::Manta::Config->create(
   version => $version,
   bams => \@bams,
   reference_fasta => $fasta,
   working_directory => $tmp_dir,
);

my $cmdline_array_ref = $config_manta->build_cmdline_array_ref();
is_deeply($cmdline_array_ref,\@expected_cmdline_array,'expected command-line array ref for version');

$config_manta->_resolve_input_files;
my @input_files = $config_manta->input_files;                       
is_deeply(\@input_files,\@expected_input_files,'expected input files');

# TODO: Test the execute subroutine
    
