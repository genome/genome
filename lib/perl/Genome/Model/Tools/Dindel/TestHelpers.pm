package Genome::Model::Tools::Dindel::TestHelpers;

use strict;
use warnings;

use Genome::Utility::Test;
use Test::More;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw(
    get_test_dir
    get_ref_fasta
    compare_output_to_test_data
);

my $SHARED_DATA_VERSION = 'v1';

sub get_test_dir {
    my ($package, $version) = @_;
    return Genome::Utility::Test->data_dir($package, "v$version");
}


sub get_ref_fasta {
    return Genome::Utility::Test->shared_test_data('human_reference_37.fa', $SHARED_DATA_VERSION);
}

sub compare_output_to_test_data {
    my ($output_path, $output_dir, $test_dir) = @_;
    Genome::Utility::Test->compare_to_blessed_file(
        output_path => $output_path,
        output_dir => $output_dir,
        test_dir => $test_dir,
    );
}
