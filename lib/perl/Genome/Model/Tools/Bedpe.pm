package Genome::Model::Tools::Bedpe;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Bedpe {
    is => ['Command::Tree'],

};

sub help_brief {
    "Tools related to BEDPE files.",
}

sub help_synopsis {
    return <<"EOS"
  gmt bedpe...
EOS
}

sub help_detail {
    return <<EOS
    This is a collection of small tools related to BEDPE files.
EOS
}

1;
