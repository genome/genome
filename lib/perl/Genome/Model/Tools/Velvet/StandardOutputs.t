#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::More;
use above "Genome";
require File::Compare;

use_ok ('Genome::Model::Tools::Velvet::StandardOutputs') or die;

my $version = 'v4';
my $data_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Velvet/StandardOutputs-'.$version;
ok (-d $data_dir, "Test suite data dir exists") or die;

#temp test dir
my $temp_dir = Genome::Sys->create_temp_directory();
ok (-d $temp_dir, "Test temp dir created") or die;

#check assembly output files exists and copy to temp dir
my @files_to_copy = qw/
    Sequences
    contigs.fa
    velvet_asm.afg
    H_GV-933124G-S.MOCK.collated.fastq
    /;
for my $file ( @files_to_copy ) {
    ok (-s $data_dir."/$file", "Test suite assembly $file file exists") or die;
    ok (File::Copy::copy( $data_dir."/$file", $temp_dir ), "Copied $file to temp test dir") or die;
}

#create/execute tool .. takes about a minute
my $create = Genome::Model::Tools::Velvet::StandardOutputs->create(
    assembly_directory => $temp_dir,
    min_contig_length => 50,
    );
ok ($create, "Created gmt velvet default tool") or die;
ok ($create->execute, "Executed gmt velvet default tool") or die;

#check output files
my @files_to_check = qw/
    contigs.bases
    contigs.quals
    gap.txt readinfo.txt
    reads.placed
    reads.unplaced
    reads.unplaced.fasta
    supercontigs.agp
    supercontigs.fasta
    /;
for my $file ( @files_to_check ) {
    ok (-e $data_dir."/edit_dir/$file", "$file exists in test dir");
    ok (-e $temp_dir."/edit_dir/$file", "$file exists in temp dir");
    ok (File::Compare::compare( $data_dir."/edit_dir/$file", $temp_dir."/edit_dir/$file") == 0, "$file files match");
}

#<STDIN>;

done_testing();
