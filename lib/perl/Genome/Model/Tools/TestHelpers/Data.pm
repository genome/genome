package Genome::Model::Tools::TestHelpers::Data;

use strict;
use warnings;

use Params::Validate qw(:types);
use above 'Genome';
use Genome::Model::Tools::TestHelpers::General qw(
    get_test_data
);
use Test::More;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw(
    get_shared_test_data
);

sub get_shared_test_data {
    my ($data_path, $version) = Params::Validate::validate_pos(@_,
        1, 1);

    return get_test_data('Genome::Model::Tools::TestHelpers::Data',
        $data_path, $version);
}

1;
