package Genome::Test::Factory::Qc::Result;

use Genome::Test::Factory::Base;
@ISA = (Genome::Test::Factory::Base);

use strict;
use warnings;
use Genome;

sub generate_obj {
    my $self = shift;
    return Genome::Qc::Result->__define__(@_);
}


1;
