#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More;
require File::Compare;

use_ok ('Genome::Model::Tools::Soap::FastaToAgp') or die;

#check test data dir and files
my $version = 'v1';
my $data_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Soap/FastaToAgp_'.$version;
ok (-d $data_dir, "Data dir exists") or die;

#check test input/output files
ok (-s $data_dir.'/TEST.scafSeq', "TEST.scafSeq file exists in data dir") or die;
my @out_files = qw/ SRS012663_PGA.scaffolds.fa SRS012663_PGA.contigs.fa SRS012663_PGA.agp /;
foreach (@out_files) {
    ok (-s $data_dir."/$_", "Test $_ file exists") or die;
}

#create tmp test dir
my $temp_dir = Genome::Sys->create_temp_directory();
ok (-d $temp_dir, "Temp test dir created") or die;

#copy input files over to temp dir
ok (File::Copy::copy ($data_dir.'/TEST.scafSeq', $temp_dir), "Copied TEST.scafSeq to temp test dir") or die;
ok (-s $temp_dir.'/TEST.scafSeq', "Temp test dir TEST.scafSeq file exists") or die;

#create/execute
my $c = Genome::Model::Tools::Soap::FastaToAgp->create (
    assembly_directory => $temp_dir,
    scaffold_size_cutoff => 100,
    output_dir => $temp_dir,
    file_prefix => 'SRS012663_PGA',
    version => '9.27.10',
    );
ok ($c, "Created fasta-to-agp tool") or die;
ok ($c->execute, "Successfully executed fasta-to-agp tool") or die;

#TODO - compare output files

done_testing();
