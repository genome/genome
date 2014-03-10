#!/gsc/bin/perl

use strict;
use warnings;

use above "Genome";
use Test::More;
use Genome::Utility::Test;

my $class = "Genome::Model::Tools::Tcga::FixSnvVcfNullAlt";
use_ok($class);

my $base_dir = Genome::Utility::Test->data_dir_ok($class, "v1");

my $input_file    = $base_dir.'/input.vcf';
my $expected_file = $base_dir.'/output.vcf';
my $output_file   = Genome::Sys->create_temp_file_path;

my $cmd = Genome::Model::Tools::Tcga::FixSnvVcfNullAlt->create(
    input_file  => $input_file,
    output_file => $output_file,
);

ok($cmd, "Command created");
ok($cmd->execute, "Command executed");

Genome::Utility::Test::compare_ok($output_file, $expected_file, "Output file was created as expected");

done_testing;
