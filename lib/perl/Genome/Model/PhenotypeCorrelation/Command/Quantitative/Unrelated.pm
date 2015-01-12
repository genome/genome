package Genome::Model::PhenotypeCorrelation::Command::Quantitative::Unrelated;

use strict;
use warnings;

use Carp 'confess';

class Genome::Model::PhenotypeCorrelation::Command::Quantitative::Unrelated {
    is  => 'Genome::Model::PhenotypeCorrelation::Command::Base',
};

sub help_synopsis {
    return <<"EOS"
EOS
}

sub help_detail {
    return <<"EOS"
Need documenation here.
EOS
}

sub execute {
    my $self = shift;

    return 1;
}

1;
