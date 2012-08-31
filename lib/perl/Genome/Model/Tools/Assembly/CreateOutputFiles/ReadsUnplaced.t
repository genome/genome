#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";
use Test::More;
require File::Compare;

use_ok ('Genome::Model::Tools::Assembly::CreateOutputFiles::ReadsUnplaced') or die;

#check data dir and input files for test
my $data_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Assembly/CreateOutputFiles/ReadsUnplaced_v0';
ok(-d $data_dir, "Data dir exists") or die;
for my $file (qw/ GABJJ9O01.fasta.gz GABJJ9O01.fasta.qual.gz GABJJ9O02.fasta.gz reads.placed /) {
    ok(-s $data_dir."/edit_dir/$file", "$file file exists in data dir");
}

#create temp directory
my $temp_dir = Genome::Sys->create_temp_directory();
Genome::Sys->create_directory( $temp_dir.'/edit_dir' );

#copy data files;
for my $file (qw/ GABJJ9O01.fasta.gz GABJJ9O01.fasta.qual.gz GABJJ9O02.fasta.gz reads.placed /) {
    ok(File::Copy::copy($data_dir."/edit_dir/$file", $temp_dir.'/edit_dir'), "Copied $file file to temp dir");
    ok(-s $temp_dir."/edit_dir/$file", "$file file exists in temp edit_dir");
}

#create / execute tool
my $c = Genome::Model::Tools::Assembly::CreateOutputFiles::ReadsUnplaced->create(
    directory => $temp_dir,
    reads_placed_file => $temp_dir.'/edit_dir/reads.placed',
    );
ok($c, "Created reads-unplaced");
ok($c->execute, "Executed reads-unplaced");

#compare output files
foreach ('reads.unplaced', 'reads.unplaced.fasta') {
    ok(-s $data_dir."/edit_dir/$_", "Data dir $_ file exists");
    ok(-s $temp_dir."/edit_dir/$_", "Temp dir $_ file exists");
    ok(File::Compare::compare($data_dir."/edit_dir/$_", $temp_dir."/edit_dir/$_") == 0, "$_ files match");
}

done_testing();

exit;
