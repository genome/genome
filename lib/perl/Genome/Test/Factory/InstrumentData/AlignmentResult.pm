package Genome::Test::Factory::InstrumentData::AlignmentResult;
use base qw(Genome::Test::Factory::Base);

use strict;
use warnings;

use Genome;

sub generate_obj {
    my $self = shift;
    return Genome::InstrumentData::AlignmentResult::Bwa->__define__(@_);
}


1;
