#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";
use Test::More;

use_ok ('Genome::Model::Tools::Assembly::ReScaffoldMsiAce');

my $data_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Assembly/ReScaffoldMsiAce_v5';
ok(-d $data_dir, "Found data dir");

my $temp_dir = Genome::Sys->create_temp_directory();
ok(-d $temp_dir, "Made temp test dir");

mkdir $temp_dir.'/edit_dir';
ok(-d $temp_dir.'/edit_dir', "Made edit_dir in temp dir");

mkdir $temp_dir.'/phdball_dir';
ok(-d $temp_dir.'/phdball_dir', "Made phdball_dir in temp dir");

ok( File::Copy::copy( $data_dir.'/phdball_dir/phd.ball.1', $temp_dir.'/phdball_dir' ), "Copied phdball file" );

foreach ('test.ace', 'scaffolds') {
    ok(-s $data_dir."/edit_dir/$_", "Data dir $_ file exists");
    ok(File::Copy::copy($data_dir."/edit_dir/$_", $temp_dir."/edit_dir/"), "Copied $_ to temp_dir");
    ok(-s $temp_dir."/edit_dir/$_", "Temp dir $_ file exists");
}

my $create = Genome::Model::Tools::Assembly::ReScaffoldMsiAce->create (
    acefile => $temp_dir.'/edit_dir/test.ace',
    scaffold_file => $temp_dir.'/edit_dir/scaffolds',
    assembly_directory => $temp_dir,
    min_contig_length => 50,
    );
ok( $create, "Created re-scaffold-msi-ace");

ok( $create->execute, "Executed re-scaffold-msi-ace successfully");

ok(-s $temp_dir.'/edit_dir/ace.msi', "Created new scaffolded ace file");
my @diff = `sdiff -s $temp_dir/edit_dir/ace.msi $data_dir/edit_dir/ace.msi`;
is(scalar (@diff), 1, "New ace file matches test ace file except for location of phdball file");

ok(-s $temp_dir.'/edit_dir/msi.gap.txt', "Created msi.gap.txt file");
my @diff2 = `sdiff -s $temp_dir/edit_dir/msi.gap.txt $data_dir/edit_dir/msi.gap.txt`;
is(scalar (@diff2), 0, "New msi.gap.txt file matches test msi.gap.txt file");

#<STDIN>;

done_testing();
