#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

require File::Compare;
require File::Temp;
use Test::More;

use_ok('Genome::Model::Tools::PooledBac::MapContigsToAssembly') or die;

my $dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-PooledBac-MapContigsToAssembly/v1';
ok(-d $dir, 'test dir exists');
my $example_contig_map = $dir.'/CONTIG_MAP';
ok(-s $example_contig_map, 'example contig map exists');

my $tmpdir = File::Temp::tempdir(CLEANUP => 1);
my $blast_output_base = 'bac_region_db.blast';
symlink("$dir/$blast_output_base", "$tmpdir/$blast_output_base");
ok(-s "$tmpdir/$blast_output_base", 'symlink blast output');

my $map_contigs = Genome::Model::Tools::PooledBac::MapContigsToAssembly->create(
    pooled_bac_dir => $dir,
    ace_file_name => 'contig_names_only.ace',
    project_dir => $tmpdir,
    #project_dir => $dir,
);
ok($map_contigs, 'create');
ok($map_contigs->execute, 'execute');

is(File::Compare::compare($example_contig_map, $tmpdir.'/CONTIG_MAP'), 0, 'contig map matches');

#print STDERR "$tmpdir\n"; <STDIN>;
done_testing();
exit;

