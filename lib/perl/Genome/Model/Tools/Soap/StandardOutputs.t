#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::More;
use above 'Genome';
require File::Compare;

use_ok ('Genome::Model::Tools::Soap::StandardOutputs') or die;

#check for test suite dir
my $data_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model/DeNovoAssembly/soap_solexa_build_v8';
ok (-d $data_dir, "Test suite data dir exists") or die;

#create temp test dir .. tool creates edit_dir above it
my $temp_dir = Genome::Sys->create_temp_directory();
ok (-d $temp_dir, "Temp test dir created") or die;

#copy post assembly input files
my @files = qw/ H_KT-185-1-0089515594_WUGC.2852968107.reverse.fastq
                H_KT-185-1-0089515594_WUGC.2852968107.forward.fastq
                H_KT-185-1-0089515594_WUGC.scafSeq /;
foreach ( @files ) {
    ok (-s $data_dir."/$_", "$_ file exists in test suite data dir") or die;
    ok (File::Copy::copy( $data_dir."/$_", $temp_dir ), "Copied $_ file to temp dir" ) or die;
}

my $tool = Genome::Model::Tools::Soap::StandardOutputs->create(
    assembly_directory => $temp_dir,
    min_contig_length => 50,
);
ok ($tool, "Created soap post assemble default tool") or die;
ok ($tool->execute, "Successfully executed soap default post assemble tool") or die;

#compare output files - these files differ with each run

done_testing();

exit;

