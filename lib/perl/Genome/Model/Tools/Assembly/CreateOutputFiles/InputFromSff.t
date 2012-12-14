#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";
use Test::More;

my $archos = `uname -a`;
unless ($archos =~ /64/) {
    plan skip_all => "Must run from 64-bit machine";
}

use_ok ('Genome::Model::Tools::Assembly::CreateOutputFiles::InputFromSff') or die;

#check test dir/files
my $data_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Assembly/CreateOutputFiles/InputFromSff_v0';
ok (-d $data_dir, "Data dir exists");
ok (-s $data_dir.'/sff/GABJJ9O01.sff', "Test sff file exists");

#create temp test dir
my $temp_dir = Genome::Sys->create_temp_directory();
my $temp_sff_dir = Genome::Sys->create_directory($temp_dir.'/sff');

#copy data over
ok (File::Copy::copy($data_dir.'/sff/GABJJ9O01.sff', $temp_dir.'/sff'), "Copied data file to temp dir");

#create/execute tool
my $create = Genome::Model::Tools::Assembly::CreateOutputFiles::InputFromSff->create(
    directory => $temp_dir,
    );
ok ($create, "Created input-from-sff tool");
ok ($create->execute, "Executed input-from-sff tool");

#compare output files
for ('GABJJ9O01.fasta.gz', 'GABJJ9O01.fasta.gz') {
    ok( -s $data_dir."/edit_dir/$_", "Data dir $_ file exists" );
    ok( -s $temp_dir."/edit_dir/$_", "New $_ file created" );
    my @diff = `zdiff $data_dir/edit_dir/$_ $temp_dir/edit_dir/$_`;
    is (scalar @diff, 0, "Data and temp $_ files match");
}

done_testing();
