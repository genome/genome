#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Test::More;
require File::Compare;

#check test suite dir/files
my $test_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Assembly-SplitScaffold';
ok (-d $test_dir, "Test suite dir exists");
foreach (qw/ merge.ace out.ace /) {
    ok( -s $test_dir."/$_", "Test suite $_ file exists");
}

#create temp test dir
my $temp_dir = Genome::Sys->create_temp_directory();
ok (-d $temp_dir, "Temp test dir created");

#create/run tool
my $create = Genome::Model::Tools::Assembly::MergeScaffolds->create(
    ace_file => $test_dir.'/merge.ace',
    left_scaffold => 'Contig60.1',
    right_scaffold => 'Contig120.1',
    out_file_name => $temp_dir.'/out.ace',
    );
ok ($create, "Successfully created merge-scaffold tool");
ok ($create->execute, "Successfully executed merge-scaffold tool");

#compare output files
ok (File::Compare::compare( $test_dir.'/out.ace', $temp_dir.'/out.ace') == 0, "Test outout files match");

done_testing();
