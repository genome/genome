package Genome::Model::Tools::Dindel::TestHelpers;

use strict;
use warnings;

use Genome::Model::Tools::TestHelpers::General qw(
    get_test_dir
    ensure_file
    compare_to_blessed_file
);
use Genome::Model::Tools::TestHelpers::Data qw(
    get_shared_test_data
);
use Test::More;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw(
    get_test_dir
    get_ref_fasta
    compare_output_to_test_data
);

my $SHARED_DATA_VERSION = 1;

sub get_ref_fasta {
    return get_shared_test_data('human_reference_37.fa', $SHARED_DATA_VERSION);
}

sub compare_output_to_test_data {
    my ($output_path, $output_dir, $test_dir) = @_;
    return compare_to_blessed_file(
        output_path => $output_path,
        output_dir => $output_dir,
        test_dir => $test_dir,
    );
}
