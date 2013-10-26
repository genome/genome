package Genome::Test::Factory::InstrumentData::Solexa;
use base qw(Genome::Test::Factory::Base);

use strict;
use warnings;

use Genome;
use Genome::Test::Factory::Library;

our @required_params = qw(library_id);

sub generate_obj {
    my $self = shift;
    return Genome::InstrumentData::Solexa->create(@_);
}

sub create_library_id {
    my $lib = Genome::Test::Factory::Library->setup_object();
    return $lib->id;
}

1;
