package Genome::Model::Tools::Bed::Convert::Snv;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Bed::Convert::Snv {
    is => ['Genome::Model::Tools::Bed::Convert'],

};

sub help_brief {
    "Tools to convert other SNV formats to BED.",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
  gmt bed convert snv...
EOS
}

sub help_detail {                           
    return <<EOS
    This is a collection of small tools to take SNV calls in various formats and convert them to a common BED format (using the first four columns).
EOS
}

1;
