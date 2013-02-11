#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Test::More tests => 4;
use File::Compare;

BEGIN {
      use_ok('Genome::Model::Tools::RepeatMasker::MergeTables');
};
my $root_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-RepeatMasker/MergeTables';
my @lanes = qw(3 4);
my @tables;
for my $lane (@lanes) {
    push @tables, $root_dir .'/s_'. $lane .'_sequence.fa.tsv';
}
my $expected_file = $root_dir .'/all.tsv';


my $tmp_dir = File::Temp::tempdir('RepeatMasker-MergeTables-'. Genome::Sys->username .'-XXXX',CLEANUP=>1, TMPDIR => 1);
my $merged_file = $tmp_dir .'/merged.tsv';
my $merge = Genome::Model::Tools::RepeatMasker::MergeTables->create(
    input_tables => \@tables,
    output_table => $merged_file,
);
isa_ok($merge,'Genome::Model::Tools::RepeatMasker::MergeTables');
ok($merge->execute,'execute command '. $merge->command_name);
ok(!compare($merged_file,$expected_file),'files are identical');
