package Genome::Model::Tools::Dindel::TestHelpers;

use strict;
use warnings;

use Genome::Utility::Test 'compare_ok';
use Test::More;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw(
    get_test_dir
    get_ref_fasta
    compare_output_to_test_data
);

sub get_test_dir {
    my ($class, $version) = @_;

    my $test_dir = File::Spec->join(Genome::Utility::Test->data_dir($class), "v$version");
    diag "Test data located at $test_dir\n";
    return $test_dir;
}

sub get_ref_fasta {
    my $reference_sequence_build = Genome::Model::Build->get(106942997);
    my $ref_fasta = $reference_sequence_build->full_consensus_path('fa');
    diag "Reference sequence fasta: $ref_fasta\n";
    ok(-s $ref_fasta, "Found Reference Sequence fasta") || die;

    return $ref_fasta;
}

sub compare_output_to_test_data {
    my ($output_path, $output_dir, $test_dir) = @_;
    (my $expected_path = $output_path) =~ s/$output_dir/$test_dir/;
    ok(-f $expected_path, "Found test-data: $expected_path") || die;
    ok(-f $output_path, "Found output-data: $output_path") || die;
    compare_ok($expected_path, $output_path, "test-data matches output-data");
}
