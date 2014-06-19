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
    test_for_version($version, $allele, $test_dir);
}

my @invalid_alleles  = qw(HLA-A02:01 A0201);
my $iterator_2 = each_array( @versions, @invalid_alleles );
while ( my ($version, $allele) = $iterator_2->() ) {
    test_for_valid_allele($version, $allele);
}

sub test_for_valid_allele {
    my ($version, $allele) = @_;

    ok(!$class->is_valid_allele_for_netmhc_version($allele, $version), "Inputs not valid");
}

sub test_for_version {
    my ($version, $allele, $test_dir) = @_;

    my $expected_output = File::Spec->join($test_dir, "expected.$version.xls");
    my $output_dir = Genome::Sys->create_temp_directory;

    my $cmd = $class->create(
        allele => $allele,
        fasta_file => $input_fasta,
        output_directory => $output_dir,
        epitope_length => 9,
        netmhc_version => $version,
        sample_name => 'test',
    );
    ok($cmd, "Created command for version $version");

    ok($cmd->execute, "Command executed for version $version");

    compare_ok($cmd->output_file, $expected_output, filters => ['^NetMHC version.*'], name => "File output as expected for version $version");
}

done_testing();
