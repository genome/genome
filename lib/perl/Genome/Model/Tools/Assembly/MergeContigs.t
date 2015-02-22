#!/usr/bin/env genome-perl

use strict;
use warnings;
use above 'Genome';
use Genome;

use Genome::Model::Tools::Assembly::MergeContigs;

#use Test::More tests => 1;
use Test::More skip_all => "Does not play nice with the test harness";

my $contigs = 'merge.ace contig00012.0 merge.ace contig00013.1';
my $path = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Assembly-MergeContigs/edit_dir';

my $output_file_name = 'out.ace';
chdir($path);
system "/bin/rm -f *.db";
ok(Genome::Model::Tools::Assembly::MergeContigs->execute(contigs => $contigs, o => $output_file_name, cc => 1)->result, "MergeContigs executed successfully");
