package Genome::Test::Factory::InstrumentData::Imported;
use base qw(Genome::Test::Factory::Base);

use strict;
use warnings;

use Genome;
use Genome::Test::Factory::Library;

our @required_params = qw(library_id sequencing_platform);

sub generate_obj {
    my $self = shift;
    return Genome::InstrumentData::Imported->create(@_);
}

sub create_library_id {
    my $lib = Genome::Test::Factory::Library->setup_object();
    return $lib->id;
}

sub create_sequencing_platform {
    return 'plink';
}

1;
