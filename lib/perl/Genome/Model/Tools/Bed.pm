package Genome::Model::Tools::Bed;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Bed {
    is => ['Command'],

};

sub help_brief {
    "Tools related to BED files (not to be confused with bedtools).",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
  gmt bed...
EOS
}

sub help_detail {                           
    return <<EOS
    This is a collection of small tools related to BED files.  (For the suite of tools called "bedtools", see `gmt bed-tools`.)
EOS
}

1;
