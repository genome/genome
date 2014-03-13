#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More;
use Genome::Utility::Test qw(compare_ok);
use List::MoreUtils qw( each_array );

my $class = 'Genome::Model::Tools::EpitopePrediction::RunNetmhc';
my $TEST_DATA_VERSION= 1;
use_ok($class);

my $test_dir = Genome::Utility::Test->data_dir_ok($class, $TEST_DATA_VERSION);
my $input_fasta = File::Spec->join($test_dir, "snvs_21mer.fa");

my @versions = qw(3.0 3.4);
my @alleles  = qw(A0201 HLA-A02:01);
my $iterator = each_array( @versions, @alleles );

while ( my ($version, $allele) = $iterator->() ) {
    test_for_version($version, $allele, $test_dir, $input_fasta);
}

sub test_for_version {
    my ($version, $allele, $test_dir, $netmhc_file) = @_;

    my $expected_output = File::Spec->join($test_dir, "expected.$version.xls");
    my $output_file = Genome::Sys->create_temp_file_path;
    my $stdout_file = Genome::Sys->create_temp_file_path;

    my $cmd = $class->create(
        allele => $allele,
        fasta_file => $input_fasta,
        output_file => $output_file,
        epitope_length => 9,
        stdout_file => $stdout_file,
        version => $version,
    );
    ok($cmd, "Created command for version $version");

    ok($cmd->execute, "Command executed for version $version");

    compare_ok($output_file, $expected_output, filters => ['^NetMHC version.*'], name => "File output as expected for version $version");
}

done_testing();
