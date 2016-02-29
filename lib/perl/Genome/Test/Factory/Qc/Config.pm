package Genome::Test::Factory::Qc::Config;

use Genome::Test::Factory::Base;
@ISA = (Genome::Test::Factory::Base);

use strict;
use warnings;
use Genome;

sub generate_obj {
    my $self = shift;
    return Genome::Qc::Config->__define__(@_);
}


1;
