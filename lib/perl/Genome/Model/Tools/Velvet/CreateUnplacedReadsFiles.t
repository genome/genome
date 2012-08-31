#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";
use Test::More;

require File::Compare;

use_ok( 'Genome::Model::Tools::Velvet::CreateUnplacedReadsFiles' );

#TODO - move to correct test suite module dir when all tests are configured
my $module = 'Genome-Model-Tools-Assembly-CreateOutputFiles2';
my $data_dir = $ENV{GENOME_TEST_INPUTS} . "/$module";

ok(-d $data_dir, "Found data directory: $data_dir");

my $temp_dir = Genome::Sys->create_temp_directory();

#link project dir files
foreach ('velvet_asm.afg', 'Sequences') {
    ok(-s $data_dir.'/'.$_, "Data dir $_ file exists"); 
    symlink ($data_dir.'/'.$_, $temp_dir.'/'.$_);
    ok(-s $temp_dir.'/'.$_, "Tmp dir $_ file exists");
}

#create/execute tool
my $create = Genome::Model::Tools::Velvet::CreateUnplacedReadsFiles->create(
    assembly_directory => $temp_dir,
    min_contig_length => 50,
    );
ok( $create, "Created tool");
ok( $create->execute, "Successfully executed tool");

foreach ('reads.unplaced', 'reads.unplaced.fasta') {
    ok(-s $data_dir."/edit_dir/$_", "Data dir $_ file exists");
    ok(-s $temp_dir."/edit_dir/$_", "Tmp dir $_ file got created");
    ok(File::Compare::compare($data_dir."/edit_dir/$_", $temp_dir."/edit_dir/$_") == 0, "$_ files match");
}

done_testing();

exit;
