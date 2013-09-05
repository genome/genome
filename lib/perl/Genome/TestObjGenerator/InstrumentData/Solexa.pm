package Genome::TestObjGenerator::InstrumentData::Solexa;
use base qw(Genome::TestObjGenerator::Base);

use strict;
use warnings;
use Genome;

sub generate_obj {
    my $self = shift;
    return Genome::InstrumentData::Solexa->create(@_);
}

1;
