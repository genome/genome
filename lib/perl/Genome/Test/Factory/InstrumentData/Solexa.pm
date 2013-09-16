package Genome::Test::Factory::InstrumentData::Solexa;
use base qw(Genome::Test::Factory::Base);

use strict;
use warnings;
use Genome;

sub generate_obj {
    my $self = shift;
    return Genome::InstrumentData::Solexa->create(@_);
}

1;
