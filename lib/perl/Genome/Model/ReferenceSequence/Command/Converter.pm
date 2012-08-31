package Genome::Model::ReferenceSequence::Command::Converter;

use strict;
use warnings;

use Genome;


class Genome::Model::ReferenceSequence::Command::Converter {
    is => 'Genome::Command::Base',
};

sub help_brief {
    "work with reference-sequence converters"
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
 gmt model reference-sequence converter ...
EOS
}

sub help_detail {
    return <<EOS
A collection of commands to interact with reference-sequence converters.
EOS
}

1;
