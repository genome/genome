package Genome::Test::Factory::Sample;
use base qw(Genome::Test::Factory::Base);

use strict;
use warnings;

use Genome;
use Genome::Test::Factory::Individual;

our @required_params = qw(name source_id);

sub generate_obj {
    my $self = shift;
    return Genome::Sample->create(@_);
}

sub create_name {
    return Genome::Test::Factory::Util::generate_name('sample_name');
}

sub create_source_id {
    my $source = Genome::Test::Factory::Individual->setup_object();
    return $source->id;
}

1;
