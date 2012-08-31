#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Test::More;
use File::Compare;

if ($] < 5.010) {
    plan skip_all => "this test is only runnable on perl 5.10+"
}
plan tests => 8;

use_ok('Genome::Model::Tools::BioSamtools::ListChromosomes');
my $data_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-BioSamtools-ListChromosomes';
my $seq_dict = $data_dir .'/seqdict.sam';

my $expected_output_file = $data_dir .'/expected_output/list.txt';

# Run once without an output file
my $list = Genome::Model::Tools::BioSamtools::ListChromosomes->create(
    input_file => $seq_dict,
);
isa_ok($list,'Genome::Model::Tools::BioSamtools::ListChromosomes');
ok($list->execute,'execute command '. $list->command_name);
is(scalar(@{$list->chromosome_array_ref}),84,'Found all 84 chromosomes');

# Run the second time printing the names to an output file
my $tmp_file = Genome::Sys->create_temp_file_path;
$list = Genome::Model::Tools::BioSamtools::ListChromosomes->create(
    input_file => $seq_dict,
    output_file => $tmp_file,
);
isa_ok($list,'Genome::Model::Tools::BioSamtools::ListChromosomes');
ok($list->execute,'execute command '. $list->command_name);
is(scalar(@{$list->chromosome_array_ref}),84,'Found all 84 chromosomes');
ok(!(compare($list->output_file,$expected_output_file)),'Output matches expected file');


exit;
