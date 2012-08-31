package Genome::Info::BamFlagstat;
use strict;
use warnings;
use Genome;

# TODO: refactor calls to this module to ::Sam::Flagstat

use Genome::Model::Tools::Sam::Flagstat;
 
sub get_data {
    my ($class, $flag_file) = @_;
    return Genome::Model::Tools::Sam::Flagstat->parse_file_into_hashref($flag_file);
}

1;
