package Genome::Model::Tools::Dindel::TestHelpers;

use strict;
use warnings;

use Genome::Model::Tools::TestHelpers::General qw(
    get_test_dir
    ensure_file
    compare_to_blessed_file
);
use Test::More;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw(
    get_test_dir
    get_ref_fasta
    compare_output_to_test_data
);

sub get_ref_fasta {
    my $reference_sequence_build = Genome::Model::Build->get(106942997);
    my $ref_fasta = $reference_sequence_build->full_consensus_path('fa');

    ensure_file($ref_fasta);
    note "Reference sequence fasta found at ($ref_fasta)\n";

    return $ref_fasta;
}

sub compare_output_to_test_data {
    my ($output_path, $output_dir, $test_dir) = @_;
    return compare_to_blessed_file(
        output_path => $output_path,
        output_dir => $output_dir,
        test_dir => $test_dir,
    );
}
